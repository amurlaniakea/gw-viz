"""Signal metrics for GW strain.

NOTE: nothing here is matched-filter SNR. The statistic is an *excess-power
significance* measured against off-source data from the same file. Published
SNRs (e.g. ~24 network SNR for GW150914) come from matched filtering and will
not match these numbers.
"""
from __future__ import annotations

import numpy as np
from scipy.signal import butter, sosfiltfilt, welch
from scipy.signal.windows import tukey

BAND = (35.0, 350.0)  # Hz, where GW150914-like mergers sit in LIGO's sensitive band


def whiten_bandpass(
    strain: np.ndarray, fs: float, band: tuple[float, float] = BAND, seglen_s: float = 4.0
) -> np.ndarray:
    """Whiten with a Welch PSD estimated from the data itself, then band-pass."""
    n = len(strain)
    # fixed 1 s taper at each end (independent of file length)
    x = (strain - np.mean(strain)) * tukey(n, alpha=min(1.0, 2.0 * fs / n))
    freqs, psd = welch(x, fs=fs, nperseg=min(int(seglen_s * fs), n), average="median")  # robust to the signal itself
    f = np.fft.rfftfreq(n, d=1.0 / fs)
    asd = np.sqrt(np.interp(f, freqs, psd))
    w = np.fft.irfft(np.fft.rfft(x) / asd, n=n)
    sos = butter(4, band, btype="bandpass", fs=fs, output="sos")
    return sosfiltfilt(sos, w)


def analyze_signal(
    times: np.ndarray,
    strain: np.ndarray,
    event_gps: float | None = None,
    window_s: float = 0.3,
    exclude_s: float = 1.0,
) -> dict:
    """On-source vs off-source excess-power statistics around `event_gps`.

    If `event_gps` is None the loudest whitened sample is used as the on-source centre.
    """
    if len(times) < 2:
        raise ValueError("Insufficient data points.")
    dt = float(times[1] - times[0])
    fs = 1.0 / dt
    w = whiten_bandpass(strain, fs)

    n_win = max(1, round(window_s * fs))
    centre = event_gps if event_gps is not None else float(times[np.argmax(np.abs(w))])
    on = (times >= centre - 0.75 * window_s) & (times < centre + 0.25 * window_s)
    if on.sum() < n_win // 2:
        raise ValueError("Event time falls outside the data span.")

    # Off-source windows, skipping 1 s at both file edges (filter transients) and around the event.
    edge = int(2 * fs)  # skip the 1 s tapers plus filter transients
    starts = np.arange(edge, len(w) - n_win - edge, n_win)
    off_e = []
    for s in starts:
        t_mid = times[s + n_win // 2]
        if abs(t_mid - centre) > exclude_s:
            off_e.append(float(np.sum(w[s : s + n_win] ** 2)))
    if len(off_e) < 5:
        raise ValueError("Not enough off-source data for a background estimate.")
    off_e = np.asarray(off_e)
    e_on = float(np.sum(w[on] ** 2))
    sigma_w = float(np.sqrt(np.mean(off_e) / n_win))

    # Peak of the band-limited (unwhitened) strain inside the window, not of raw strain.
    sos = butter(4, BAND, btype="bandpass", fs=fs, output="sos")
    bp = sosfiltfilt(sos, strain - np.mean(strain))
    peak_amp = float(np.max(np.abs(bp[on])))

    # Power-weighted central frequency inside the window.
    seg = bp[on] * np.hanning(on.sum())
    spec = np.abs(np.fft.rfft(seg)) ** 2
    ff = np.fft.rfftfreq(len(seg), d=dt)
    m = (ff >= BAND[0]) & (ff <= BAND[1])
    central = float(np.sum(ff[m] * spec[m]) / np.sum(spec[m])) if np.any(spec[m]) else 0.0

    return {
        "excess_power_sigma": float((e_on - off_e.mean()) / off_e.std(ddof=1)),
        "peak_whitened_sigma": float(np.max(np.abs(w[on])) / sigma_w),
        "peak_amp_bandpassed": peak_amp,
        "central_freq_hz": central,
        "n_off_windows": len(off_e),
        "data_span_s": float(times[-1] - times[0]),
        "sample_rate_hz": float(fs),
        "event_gps": float(centre),
    }
