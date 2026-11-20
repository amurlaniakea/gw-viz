/* GW Analysis Toolkit - demo didáctica en el navegador.
 *
 * DATOS SINTÉTICOS: ruido gaussiano coloreado + una señal de inspiral inyectada. Reproduce la misma
 * matemática que src/gw_analysis/processing/matched_filter.py (plantillas TaylorF2 2PN cortadas en el
 * ISCO, PSD de Welch con mediana, SNR maximizado en fase, fondo en ventanas sin señal).
 * No descarga ni analiza datos reales: eso lo hace el paquete Python (ver README).
 * Licencia: AGPL-3.0-or-later. (c) Pedro Sordo Martínez.
 */
(function (root) {
  'use strict';

  var MSUN_S = 4.925491025543576e-6; // G * Msun / c^3 (s)
  var FS = 1024;                     // Hz
  var DUR = 32;                      // s
  var N = FS * DUR;                  // 32768 (potencia de 2)
  var EVENT_T = DUR / 2;             // instante de coalescencia de la plantilla inyectada (s)

  /* ---------- utilidades numéricas ---------- */
  function mulberry32(seed) {
    var a = seed >>> 0;
    return function () {
      a = (a + 0x6D2B79F5) | 0;
      var t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  function gauss(rand) {
    var u = 0;
    while (u === 0) u = rand();
    var v = rand();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }

  function fftInPlace(re, im) {
    var n = re.length, i, j, bit, t;
    for (i = 1, j = 0; i < n; i++) {
      bit = n >> 1;
      for (; j & bit; bit >>= 1) j ^= bit;
      j ^= bit;
      if (i < j) {
        t = re[i]; re[i] = re[j]; re[j] = t;
        t = im[i]; im[i] = im[j]; im[j] = t;
      }
    }
    for (var len = 2; len <= n; len <<= 1) {
      var ang = -2 * Math.PI / len, wr = Math.cos(ang), wi = Math.sin(ang), half = len >> 1;
      for (i = 0; i < n; i += len) {
        var cr = 1, ci = 0;
        for (var k = 0; k < half; k++) {
          var a = i + k, b = a + half;
          var vr = re[b] * cr - im[b] * ci, vi = re[b] * ci + im[b] * cr;
          re[b] = re[a] - vr; im[b] = im[a] - vi;
          re[a] += vr; im[a] += vi;
          var nr = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = nr;
        }
      }
    }
  }

  function ifftInPlace(re, im) { // escalada 1/n, como numpy.fft.ifft
    var n = re.length, i;
    for (i = 0; i < n; i++) im[i] = -im[i];
    fftInPlace(re, im);
    for (i = 0; i < n; i++) { re[i] /= n; im[i] = -im[i] / n; }
  }

  function rfft(x) {
    var n = x.length, re = new Float64Array(n), im = new Float64Array(n), nf = n / 2 + 1;
    re.set(x);
    fftInPlace(re, im);
    return { re: re.slice(0, nf), im: im.slice(0, nf) };
  }

  function irfft(Xre, Xim, n) {
    var re = new Float64Array(n), im = new Float64Array(n), nf = n / 2 + 1, k;
    for (k = 0; k < nf; k++) { re[k] = Xre[k]; im[k] = Xim[k]; }
    for (k = nf; k < n; k++) { re[k] = Xre[n - k]; im[k] = -Xim[n - k]; }
    im[0] = 0; im[n / 2] = 0;
    ifftInPlace(re, im);
    return re;
  }

  function tukey(m, alpha) { // scipy.signal.windows.tukey (simétrica)
    var w = new Float64Array(m), i;
    if (alpha <= 0) { w.fill(1); return w; }
    if (alpha >= 1) {
      for (i = 0; i < m; i++) w[i] = 0.5 - 0.5 * Math.cos(2 * Math.PI * i / (m - 1));
      return w;
    }
    var width = Math.floor(alpha * (m - 1) / 2);
    for (i = 0; i < m; i++) {
      if (i <= width) w[i] = 0.5 * (1 + Math.cos(Math.PI * (-1 + 2 * i / alpha / (m - 1))));
      else if (i < m - width - 1) w[i] = 1;
      else w[i] = 0.5 * (1 + Math.cos(Math.PI * (-2 / alpha + 1 + 2 * i / alpha / (m - 1))));
    }
    return w;
  }

  function median(arr) {
    var a = Array.prototype.slice.call(arr).sort(function (p, q) { return p - q; });
    var h = a.length >> 1;
    return a.length % 2 ? a[h] : 0.5 * (a[h - 1] + a[h]);
  }

  function medianBias(n) { // scipy.signal.spectral._median_bias
    var s = 1;
    for (var i = 1; i <= Math.floor((n - 1) / 2); i++) s += 1 / (2 * i + 1) - 1 / (2 * i);
    return s;
  }

  /* PSD unilateral de Welch con promedio por mediana (scipy.signal.welch, average='median') */
  function welchMedian(x, fs, nperseg) {
    var n = x.length, step = nperseg >> 1, nseg = Math.floor((n - nperseg) / step) + 1;
    var win = new Float64Array(nperseg), s2 = 0, i, k, s;
    for (i = 0; i < nperseg; i++) { win[i] = 0.5 - 0.5 * Math.cos(2 * Math.PI * i / nperseg); s2 += win[i] * win[i]; }
    var scale = 1 / (fs * s2), nf = nperseg / 2 + 1;
    var vals = [];
    for (k = 0; k < nf; k++) vals.push(new Float64Array(nseg));
    for (s = 0; s < nseg; s++) {
      var seg = new Float64Array(nperseg), mean = 0;
      for (i = 0; i < nperseg; i++) mean += x[s * step + i];
      mean /= nperseg;
      for (i = 0; i < nperseg; i++) seg[i] = (x[s * step + i] - mean) * win[i];
      var F = rfft(seg);
      for (k = 0; k < nf; k++) {
        var p = (F.re[k] * F.re[k] + F.im[k] * F.im[k]) * scale;
        if (k > 0 && k < nf - 1) p *= 2;
        vals[k][s] = p;
      }
    }
    var bias = medianBias(nseg), freqs = new Float64Array(nf), psd = new Float64Array(nf);
    for (k = 0; k < nf; k++) { freqs[k] = k * fs / nperseg; psd[k] = median(vals[k]) / bias; }
    return { freqs: freqs, psd: psd };
  }

  function interpUniform(f, f0, df, ys) { // equivalente a numpy.interp en rejilla uniforme
    var pos = (f - f0) / df;
    if (pos <= 0) return ys[0];
    if (pos >= ys.length - 1) return ys[ys.length - 1];
    var i = Math.floor(pos), fr = pos - i;
    return ys[i] * (1 - fr) + ys[i + 1] * fr;
  }

  /* ---------- física ---------- */
  function chirpToMasses(mc, eta) {
    var total = mc * Math.pow(eta, -0.6), root = Math.sqrt(Math.max(1 - 4 * eta, 0));
    return [total * (1 + root) / 2, total * (1 - root) / 2];
  }

  function defaultBank() {
    var bank = [], etas = [0.18, 0.21, 0.23, 0.25];
    for (var mc = 15; mc <= 45.0001; mc += 2.5) {
      for (var j = 0; j < etas.length; j++) bank.push(chirpToMasses(mc, etas[j]));
    }
    return bank;
  }

  function inspiralDuration(mc, fLow) {
    var mcS = mc * MSUN_S;
    return 5 / 256 * Math.pow(mcS, -5 / 3) * Math.pow(Math.PI * fLow, -8 / 3);
  }

  /* Plantilla TaylorF2 con fase hasta 2PN, convención e^{-2 pi i f t} (NumPy). Cero fuera de [fLow, f_isco*factor]. */
  function taylorF2(freqs, m1, m2, fLow, fCutFactor) {
    var mtot = (m1 + m2) * MSUN_S, eta = m1 * m2 / Math.pow(m1 + m2, 2);
    var fCut = (fCutFactor === undefined ? 1 : fCutFactor) / (Math.pow(6, 1.5) * Math.PI * mtot);
    var n = freqs.length, re = new Float64Array(n), im = new Float64Array(n);
    for (var k = 0; k < n; k++) {
      var f = freqs[k];
      if (f < fLow || f > fCut) continue;
      var v = Math.cbrt(Math.PI * mtot * f);
      var psi = 3 / (128 * eta * Math.pow(v, 5)) * (1 +
        20 / 9 * (743 / 336 + 11 * eta / 4) * v * v -
        16 * Math.PI * v * v * v +
        10 * (3058673 / 1016064 + 5429 * eta / 1008 + 617 * eta * eta / 144) * Math.pow(v, 4));
      var amp = Math.pow(f, -7 / 6), ph = -(psi - Math.PI / 4);
      re[k] = amp * Math.cos(ph); im[k] = amp * Math.sin(ph);
    }
    return { re: re, im: im };
  }

  /* PSD analítica del ruido sintético (la misma que tests/conftest.py::colored_psd) */
  function truePsd(f, fs) { return (2 / fs) * 1e-42 / (1 + Math.pow(30 / Math.max(f, 1), 4)); }

  function coloredNoise(n, fs, seed) {
    var rand = mulberry32(seed), w = new Float64Array(n), i;
    for (i = 0; i < n; i++) w[i] = gauss(rand);
    var X = rfft(w), nf = n / 2 + 1;
    for (var k = 0; k < nf; k++) {
      var g = 1 / Math.sqrt(1 + Math.pow(30 / Math.max(k * fs / n, 1), 4));
      X.re[k] *= g; X.im[k] *= g;
    }
    var x = irfft(X.re, X.im, n);
    for (i = 0; i < n; i++) x[i] *= 1e-21;
    return x;
  }

  /* Inyecta una plantilla con SNR óptimo conocido (PSD analítica) que termina en tc segundos. */
  function injectSignal(noise, fs, tc, mc, eta, rhoOpt, fLow) {
    var n = noise.length, nf = n / 2 + 1, df = fs / n, freqs = new Float64Array(nf), k;
    for (k = 0; k < nf; k++) freqs[k] = k * df;
    var out = new Float64Array(noise), sig = new Float64Array(n);
    if (!(rhoOpt > 0)) return { data: out, signal: sig };
    var m = chirpToMasses(mc, eta), h = taylorF2(freqs, m[0], m[1], fLow === undefined ? 20 : fLow);
    var s2 = 0;
    for (k = 0; k < nf; k++) s2 += (h.re[k] * h.re[k] + h.im[k] * h.im[k]) / truePsd(freqs[k], fs);
    var sigma = Math.sqrt(4 * df * s2);
    if (!(sigma > 0)) return { data: out, signal: sig };
    var amp = rhoOpt / sigma, Sre = new Float64Array(nf), Sim = new Float64Array(nf);
    for (k = 0; k < nf; k++) {
      var ph = -2 * Math.PI * freqs[k] * tc, c = Math.cos(ph), s = Math.sin(ph);
      Sre[k] = amp * fs * (h.re[k] * c - h.im[k] * s);
      Sim[k] = amp * fs * (h.re[k] * s + h.im[k] * c);
    }
    sig = irfft(Sre, Sim, n);
    for (var i = 0; i < n; i++) out[i] += sig[i];
    return { data: out, signal: sig };
  }

  /* Filtro adaptado: mismo algoritmo que matched_filter_search() de Python. */
  function matchedFilter(x, fs, eventT, options) {
    var o = options || {};
    var fLow = o.fLow === undefined ? 20 : o.fLow, windowS = o.windowS || 0.4, excludeS = o.excludeS || 1.0;
    var bank = o.bank || defaultBank();
    var n = x.length, dt = 1 / fs, nf = n / 2 + 1, df = 1 / (n * dt), i, k;

    var mean = 0;
    for (i = 0; i < n; i++) mean += x[i];
    mean /= n;
    var taper = tukey(n, Math.min(1, 2 * fs / n)), xw = new Float64Array(n);
    for (i = 0; i < n; i++) xw[i] = (x[i] - mean) * taper[i];

    var W = welchMedian(xw, fs, Math.min(Math.floor(4 * fs), n));
    var dfw = W.freqs[1] - W.freqs[0], S = new Float64Array(nf), freqs = new Float64Array(nf);
    for (k = 0; k < nf; k++) { freqs[k] = k * df; S[k] = interpUniform(freqs[k], 0, dfw, W.psd); }
    var D = rfft(xw);
    for (k = 0; k < nf; k++) { D.re[k] *= dt; D.im[k] *= dt; }

    var rhoMax = new Float64Array(n), best = new Int32Array(n);
    var sre = new Float64Array(n), sim = new Float64Array(n);
    for (var t = 0; t < bank.length; t++) {
      var h = taylorF2(freqs, bank[t][0], bank[t][1], fLow, o.fCutFactor);
      var s2 = 0, any = false;
      for (k = 0; k < nf; k++) {
        var a2 = h.re[k] * h.re[k] + h.im[k] * h.im[k];
        if (a2 > 0) { any = true; s2 += a2 / S[k]; }
      }
      if (!any) continue;
      var sigma = Math.sqrt(4 * df * s2);
      sre.fill(0); sim.fill(0);
      for (k = 0; k < nf; k++) {
        if (h.re[k] === 0 && h.im[k] === 0) continue;
        sre[k] = (D.re[k] * h.re[k] + D.im[k] * h.im[k]) / S[k];
        sim[k] = (D.im[k] * h.re[k] - D.re[k] * h.im[k]) / S[k];
      }
      ifftInPlace(sre, sim);
      var scale = 4 * df * n / sigma;
      for (i = 0; i < n; i++) {
        var rho = scale * Math.sqrt(sre[i] * sre[i] + sim[i] * sim[i]);
        if (rho > rhoMax[i]) { rhoMax[i] = rho; best[i] = t; }
      }
    }

    var half = windowS / 2, jOn = -1;
    for (i = 0; i < n; i++) {
      if (Math.abs(i * dt - eventT) <= half && (jOn < 0 || rhoMax[i] > rhoMax[jOn])) jOn = i;
    }
    if (jOn < 0) throw new Error('El instante del evento queda fuera de los datos.');

    var nWin = Math.max(1, Math.round(windowS * fs)), edge = Math.round(2 * fs), tauMax = 0;
    for (t = 0; t < bank.length; t++) {
      var mc = Math.pow(bank[t][0] * bank[t][1], 0.6) / Math.pow(bank[t][0] + bank[t][1], 0.2);
      tauMax = Math.max(tauMax, inspiralDuration(mc, fLow));
    }
    var skipBefore = Math.max(excludeS, windowS), skipAfter = tauMax + windowS + 0.5, bg = [];
    for (var s = edge; s < n - nWin - edge; s += nWin) {
      var start = s * dt - eventT, end = (s + nWin - 1) * dt - eventT;
      if (end > -skipBefore && start < skipAfter) continue;
      var mx = 0;
      for (i = s; i < s + nWin; i++) if (rhoMax[i] > mx) mx = rhoMax[i];
      bg.push(mx);
    }
    if (bg.length < 5) throw new Error('No hay suficientes datos fuera de la fuente.');
    var ge = 0;
    for (i = 0; i < bg.length; i++) if (bg[i] >= rhoMax[jOn]) ge++;
    var m1 = bank[best[jOn]][0], m2 = bank[best[jOn]][1];

    // datos blanqueados (para dibujar): banda 20-400 Hz
    var Wre = new Float64Array(nf), Wim = new Float64Array(nf);
    for (k = 0; k < nf; k++) {
      if (freqs[k] >= 20 && freqs[k] <= 400) { var r = 1 / Math.sqrt(S[k]); Wre[k] = D.re[k] * r; Wim[k] = D.im[k] * r; }
    }
    return {
      snr: rhoMax[jOn], peakOffsetS: jOn * dt - eventT, m1: m1, m2: m2,
      chirpMass: Math.pow(m1 * m2, 0.6) / Math.pow(m1 + m2, 0.2),
      backgroundMax: Math.max.apply(null, bg), backgroundMedian: median(bg),
      falseAlarmBound: (1 + ge) / (1 + bg.length), nOffWindows: bg.length, nTemplates: bank.length,
      rho: rhoMax, whitened: irfft(Wre, Wim, n)
    };
  }

  var core = {
    FS: FS, DUR: DUR, N: N, EVENT_T: EVENT_T,
    fft: fftInPlace, ifft: ifftInPlace, rfft: rfft, irfft: irfft, tukey: tukey,
    welchMedian: welchMedian, chirpToMasses: chirpToMasses, defaultBank: defaultBank,
    taylorF2: taylorF2, truePsd: truePsd, coloredNoise: coloredNoise,
    injectSignal: injectSignal, matchedFilter: matchedFilter, mulberry32: mulberry32
  };

  /* ---------- interfaz ---------- */
  function initUI(doc, win) {
    var el = function (id) { return doc.getElementById(id); };
    var ctl = { mc: el('mc'), eta: el('eta'), snr: el('snr'), seed: el('seed') };
    var showSig = el('showSignal'), status = el('status'), stats = el('stats');
    var canv = [el('c1'), el('c2'), el('c3')];
    var cache = { seed: null, noise: null }, timer = null;

    function cssVar(name, fallback) {
      var v = '';
      try { v = win.getComputedStyle(doc.documentElement).getPropertyValue(name); } catch (e) { v = ''; }
      return (v && v.trim()) || fallback;
    }

    function plot(canvas, series, opt) {
      var dpr = win.devicePixelRatio || 1, w = canvas.clientWidth || 640, h = opt.height || 200;
      canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr);
      var ctx = canvas.getContext('2d');
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      var fg = cssVar('--fg', '#222'), grid = cssVar('--grid', '#8884'), bg = cssVar('--panel', '#fff');
      ctx.fillStyle = bg; ctx.fillRect(0, 0, w, h);
      var L = 58, R = 12, T = 22, B = 34, pw = w - L - R, ph = h - T - B;
      var x0 = opt.xmin, x1 = opt.xmax, y0 = Infinity, y1 = -Infinity, s, i;
      for (s = 0; s < series.length; s++) {
        for (i = 0; i < series[s].x.length; i++) {
          var xx = series[s].x[i];
          if (xx < x0 || xx > x1) continue;
          var yy = series[s].y[i];
          if (yy < y0) y0 = yy; if (yy > y1) y1 = yy;
        }
      }
      (opt.hlines || []).forEach(function (l) { if (l.y > y1) y1 = l.y; if (l.y < y0) y0 = l.y; });
      if (opt.zeroBased) y0 = 0;
      if (!isFinite(y0) || y0 === y1) { y0 = -1; y1 = 1; }
      var pad = (y1 - y0) * 0.08; y1 += pad; if (!opt.zeroBased) y0 -= pad;
      var X = function (v) { return L + (v - x0) / (x1 - x0) * pw; };
      var Y = function (v) { return T + (1 - (v - y0) / (y1 - y0)) * ph; };
      ctx.strokeStyle = grid; ctx.lineWidth = 1; ctx.fillStyle = fg; ctx.font = '11px sans-serif';
      ctx.textAlign = 'center';
      for (i = 0; i <= 4; i++) {
        var xv = x0 + (x1 - x0) * i / 4, px = X(xv);
        ctx.beginPath(); ctx.moveTo(px, T); ctx.lineTo(px, T + ph); ctx.stroke();
        ctx.fillText(xv.toFixed(1), px, T + ph + 14);
      }
      ctx.textAlign = 'right';
      for (i = 0; i <= 4; i++) {
        var yv = y0 + (y1 - y0) * i / 4, py = Y(yv);
        ctx.beginPath(); ctx.moveTo(L, py); ctx.lineTo(L + pw, py); ctx.stroke();
        ctx.fillText(opt.fmt ? opt.fmt(yv) : yv.toFixed(1), L - 6, py + 3);
      }
      ctx.textAlign = 'left'; ctx.fillText(opt.title, L, 13);
      ctx.textAlign = 'center'; ctx.fillText(opt.xlabel, L + pw / 2, h - 4);
      ctx.save(); ctx.beginPath(); ctx.rect(L, T, pw, ph); ctx.clip();
      series.forEach(function (sr) {
        ctx.strokeStyle = sr.color; ctx.lineWidth = sr.width || 1; ctx.beginPath();
        var started = false, col = -1, mn = 0, mx = 0, flush = function () {
          if (col < 0) return;
          var px = L + col + 0.5;
          if (!started) { ctx.moveTo(px, Y(mn)); started = true; } else ctx.lineTo(px, Y(mn));
          if (mx !== mn) ctx.lineTo(px, Y(mx));
        };
        for (var j = 0; j < sr.x.length; j++) {
          if (sr.x[j] < x0 || sr.x[j] > x1) continue;
          var c = Math.floor(X(sr.x[j]) - L), v = sr.y[j];
          if (c !== col) { flush(); col = c; mn = v; mx = v; } else { if (v < mn) mn = v; if (v > mx) mx = v; }
        }
        flush(); ctx.stroke();
      });
      (opt.hlines || []).forEach(function (l) {
        ctx.strokeStyle = l.color; ctx.setLineDash([5, 4]); ctx.beginPath();
        ctx.moveTo(L, Y(l.y)); ctx.lineTo(L + pw, Y(l.y)); ctx.stroke(); ctx.setLineDash([]);
      });
      ctx.restore();
      (opt.hlines || []).forEach(function (l) {
        ctx.fillStyle = l.color; ctx.textAlign = 'right'; ctx.fillText(l.label, L + pw - 4, Y(l.y) - 4);
      });
      if (opt.marker) {
        ctx.fillStyle = opt.marker.color; ctx.beginPath();
        ctx.arc(X(opt.marker.x), Y(opt.marker.y), 4, 0, 2 * Math.PI); ctx.fill();
      }
    }

    function setText(node, text) { if (node) node.textContent = text; }

    function compute() {
      var t0 = Date.now();
      var mc = +ctl.mc.value, eta = +ctl.eta.value, rho = +ctl.snr.value, seed = Math.floor(+ctl.seed.value) || 1;
      el('mcVal').textContent = mc.toFixed(1);
      el('etaVal').textContent = eta.toFixed(2);
      el('snrVal').textContent = rho.toFixed(0);
      if (cache.seed !== seed) { cache.noise = coloredNoise(N, FS, seed); cache.seed = seed; }
      var inj = injectSignal(cache.noise, FS, EVENT_T, mc, eta, rho, 20);
      var r = matchedFilter(inj.data, FS, EVENT_T);
      var t = new Float64Array(N), i, sd = 0, sw = 0;
      for (i = 0; i < N; i++) { t[i] = i / FS - EVENT_T; }
      for (i = 0; i < N; i++) { sd += inj.data[i] * inj.data[i]; sw += r.whitened[i] * r.whitened[i]; }
      sd = Math.sqrt(sd / N); sw = Math.sqrt(sw / N);
      var raw = new Float64Array(N), wh = new Float64Array(N), sg = new Float64Array(N);
      for (i = 0; i < N; i++) { raw[i] = inj.data[i] / 1e-21; wh[i] = r.whitened[i] / sw; sg[i] = inj.signal[i] / 1e-21; }
      var show = showSig && showSig.checked, accent = cssVar('--accent', '#0b7285'), warn = cssVar('--warn', '#c92a2a');
      var s1 = [{ x: t, y: raw, color: accent }], s2 = [{ x: t, y: wh, color: accent }];
      if (show && rho > 0) s1.push({ x: t, y: sg, color: warn, width: 1.5 });
      plot(canv[0], s1, { title: 'Datos (ruido + señal inyectada), strain x 1e-21', xlabel: 't - t_evento (s)',
        xmin: -1.5, xmax: 0.5, height: 190, fmt: function (v) { return v.toFixed(1); } });
      plot(canv[1], s2, { title: 'Datos blanqueados (20-400 Hz), unidades de desviación típica', xlabel: 't - t_evento (s)',
        xmin: -1.5, xmax: 0.5, height: 190, fmt: function (v) { return v.toFixed(1); } });
      plot(canv[2], [{ x: t, y: r.rho, color: accent }], {
        title: 'SNR del filtro adaptado (máximo sobre ' + r.nTemplates + ' plantillas)', xlabel: 't - t_evento (s)',
        xmin: -4, xmax: 4, height: 220, zeroBased: true,
        hlines: [{ y: r.backgroundMax, color: warn, label: 'máx. fondo ' + r.backgroundMax.toFixed(1) }],
        marker: { x: r.peakOffsetS, y: r.snr, color: warn }
      });
      var lines = [];
      lines.push('SNR óptimo inyectado: ' + (rho > 0 ? rho.toFixed(1) : '0 (solo ruido)'));
      lines.push('SNR recuperado: ' + r.snr.toFixed(1) + '   (mayor SNR en ' + r.nOffWindows + ' ventanas sin señal: ' + r.backgroundMax.toFixed(1) + ')');
      lines.push('Mejor plantilla: m1 = ' + r.m1.toFixed(0) + ', m2 = ' + r.m2.toFixed(0) + ' M☉ (masa chirp ' + r.chirpMass.toFixed(1) + '; inyectada ' + mc.toFixed(1) + ')');
      lines.push('Cota de falsa alarma: <= ' + r.falseAlarmBound.toFixed(3) + ' (mínimo posible con este fondo: ' + (1 / (1 + r.nOffWindows)).toFixed(3) + ')');
      setText(stats, lines.join('\n'));
      setText(status, 'Listo (' + (Date.now() - t0) + ' ms). Datos sintéticos.');
    }

    function schedule() {
      setText(status, 'Calculando...');
      if (timer) win.clearTimeout(timer);
      timer = win.setTimeout(function () {
        try { compute(); } catch (e) { setText(status, 'Error: ' + e.message); }
      }, 60);
    }

    ['mc', 'eta', 'snr', 'seed'].forEach(function (k) {
      ctl[k].addEventListener('input', schedule);
      ctl[k].addEventListener('change', schedule);
    });
    if (showSig) showSig.addEventListener('change', schedule);
    var btn = el('btnNew');
    if (btn) btn.addEventListener('click', function () {
      ctl.seed.value = String(1 + Math.floor(Math.random() * 99999)); schedule();
    });
    // figuras reales: si aún no se han generado, mostrar un aviso en vez de una imagen rota
    var imgs = doc.querySelectorAll ? doc.querySelectorAll('img[data-real-figure]') : [];
    Array.prototype.forEach.call(imgs, function (img) {
      var miss = function () {
        var p = doc.createElement('p');
        p.className = 'missing';
        p.textContent = 'Figura no disponible todavía: genérala con el comando indicado más abajo.';
        if (img.parentNode) img.parentNode.replaceChild(p, img);
      };
      img.addEventListener('error', miss);
      if (img.complete && img.naturalWidth === 0) miss();
    });
    schedule();
    return { compute: compute };
  }

  core.initUI = initUI;
  if (typeof module !== 'undefined' && module.exports) module.exports = core;
  if (typeof document !== 'undefined' && typeof window !== 'undefined') {
    root.GWDemo = core;
    var start = function () { initUI(document, window); };
    if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', start); else start();
  }
})(typeof window !== 'undefined' ? window : globalThis);
