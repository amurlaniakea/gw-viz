/* Pruebas de la demo del navegador (docs/demo.js). Ejecutar: node tests/demo.test.js */
'use strict';
const assert = require('assert');
const D = require('../docs/demo.js');

function close(a, b, tol, msg) { assert.ok(Math.abs(a - b) <= tol, `${msg}: ${a} vs ${b} (tol ${tol})`); }
let passed = 0;
function test(name, fn) { fn(); passed++; console.log('ok -', name); }

test('FFT: impulso, ida y vuelta', () => {
  const n = 64, re = new Float64Array(n), im = new Float64Array(n);
  re[0] = 1; D.fft(re, im);
  for (let k = 0; k < n; k++) { close(re[k], 1, 1e-12, 're'); close(im[k], 0, 1e-12, 'im'); }
  const x = Float64Array.from({ length: 256 }, (_, i) => Math.sin(i * 0.3) + 0.1 * i % 3);
  const X = D.rfft(x), y = D.irfft(X.re, X.im, 256);
  for (let i = 0; i < 256; i++) close(y[i], x[i], 1e-10, 'roundtrip');
});

test('FFT: una senoide cae en su bin', () => {
  const n = 1024, x = Float64Array.from({ length: n }, (_, i) => Math.cos(2 * Math.PI * 50 * i / n));
  const X = D.rfft(x);
  close(Math.hypot(X.re[50], X.im[50]), n / 2, 1e-8, 'amplitud');
  close(Math.hypot(X.re[51], X.im[51]), 0, 1e-8, 'vecino');
});

// Valores de referencia calculados con src/gw_analysis/processing/matched_filter.py (Python)
test('plantilla TaylorF2 coincide con la de Python', () => {
  const [m1, m2] = D.chirpToMasses(30.0, 0.23);
  close(m1, 46.47591715124449, 1e-9, 'm1');
  close(m2, 25.9817843258008, 1e-9, 'm2');
  const h = D.taylorF2(Float64Array.from([25, 50, 80]), m1, m2, 20);
  close(h.re[0], -0.014553553943144467, 1e-9, 're 25');
  close(h.im[0], -0.018313557015515222, 1e-9, 'im 25');
  close(h.re[1], 0.010401639684913424, 1e-9, 're 50');
  close(h.im[1], -0.000618543878015426, 1e-9, 'im 50');
  assert.strictEqual(h.re[2], 0); assert.strictEqual(h.im[2], 0); // por encima del ISCO
});

test('PSD de Welch (mediana) reproduce la PSD analítica en promedio', () => {
  // un solo bin es ruidoso (mediana de ~15 segmentos): se promedia la razón en una banda
  const x = D.coloredNoise(D.N, D.FS, 5);
  const W = D.welchMedian(x, D.FS, 4 * D.FS);
  for (const [lo, hi] of [[60, 120], [120, 250], [250, 480]]) {
    let sum = 0, cnt = 0;
    for (let k = 0; k < W.freqs.length; k++) {
      if (W.freqs[k] >= lo && W.freqs[k] < hi) { sum += W.psd[k] / D.truePsd(W.freqs[k], D.FS); cnt++; }
    }
    const ratio = sum / cnt;
    assert.ok(ratio > 0.9 && ratio < 1.1, `PSD ${lo}-${hi} Hz: razón media ${ratio}`);
  }
});

function search(rho, seed, mc, eta) {
  const noise = D.coloredNoise(D.N, D.FS, seed);
  const inj = D.injectSignal(noise, D.FS, D.EVENT_T, mc, eta, rho, 20);
  return D.matchedFilter(inj.data, D.FS, D.EVENT_T);
}

test('inyección y recuperación (SNR 25, Mc 30)', () => {
  for (const seed of [11, 21, 22]) {
    const r = search(25, seed, 30, 0.23);
    assert.ok(r.snr > 18 && r.snr < 31, `SNR recuperado ${r.snr} (seed ${seed})`);
    assert.ok(Math.abs(r.chirpMass - 30) <= 2.6, `Mc ${r.chirpMass}`);
    assert.ok(Math.abs(r.peakOffsetS) < 0.05, `offset ${r.peakOffsetS}`);
    assert.ok(Math.abs(r.falseAlarmBound - 1 / (1 + r.nOffWindows)) < 1e-12, 'FAP');
  }
});

test('solo ruido: sin picos fuertes', () => {
  for (const seed of [1, 2, 3]) {
    const r = search(0, seed, 30, 0.23);
    assert.ok(r.snr < 7.5, `SNR de ruido ${r.snr}`);
    assert.ok(r.falseAlarmBound > 0.01);
  }
});

test('una señal fuerte no contamina el fondo', () => {
  const r = search(40, 7, 25, 0.21);
  assert.ok(r.snr > 30 && r.backgroundMax < 8, `snr ${r.snr}, fondo ${r.backgroundMax}`);
});

// Prueba de humo de la interfaz con un DOM de pega: no debe lanzar excepciones y debe dibujar.
test('interfaz: se calcula y se dibuja sin errores', () => {
  const draws = { n: 0 };
  const ctx = new Proxy({}, { get: (_, p) => (p === 'canvas' ? {} : (...a) => { draws.n++; }), set: () => true });
  const nodes = {};
  const mk = (id, extra = {}) => (nodes[id] = Object.assign({
    id, value: '', textContent: '', checked: true, clientWidth: 600, width: 0, height: 0,
    addEventListener() {}, getContext: () => ctx
  }, extra));
  ['mc', 'eta', 'snr', 'seed'].forEach((k) => mk(k));
  nodes.mc.value = '30'; nodes.eta.value = '0.23'; nodes.snr.value = '25'; nodes.seed.value = '11';
  ['mcVal', 'etaVal', 'snrVal', 'status', 'stats', 'showSignal', 'btnNew', 'c1', 'c2', 'c3'].forEach((k) => mk(k));
  const doc = { getElementById: (id) => nodes[id], documentElement: {}, querySelectorAll: () => [], createElement: () => ({}) };
  const win = { devicePixelRatio: 1, getComputedStyle: () => ({ getPropertyValue: () => '' }), setTimeout: (f) => 0, clearTimeout() {} };
  const ui = D.initUI(doc, win);
  ui.compute();
  assert.ok(draws.n > 100, 'se dibujó algo');
  assert.ok(/SNR recuperado: \d/.test(nodes.stats.textContent), nodes.stats.textContent);
  assert.ok(/Listo/.test(nodes.status.textContent));
});

console.log(`\n${passed} pruebas JS correctas`);
