"""Educational matched filter with a 2PN TaylorF2 inspiral template bank.

What this is: a frequency-domain matched filter (the same idea as the LIGO/Virgo
search pipelines, e.g. FINDCHIRP) built only with NumPy/SciPy, so students can read it.

What this is NOT: parameter estimation, nor a full inspiral-merger-ringdown search.
Templates stop at the innermost stable circular orbit (ISCO), so the recovered SNR is a
*lower bound* of what a full waveform model would give, and best-fit masses are only
indicative. Published SNRs (e.g. ~24 network SNR for GW150914) will not be reproduced.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.signal import welch
from scipy.signal.windows import tukey

MSUN_S = 4.925491025543576e-6  # G * Msun / c^3, in seconds
MAX_DURATION_S = 300.0  # keep memory/time reasonable: crop the data first


def chirp_to_masses(chirp_mass: float, eta: float) -> tuple[float, float]:
    """(chirp mass, symmetric mass ratio) -> component masses in solar masses."""
    total = chirp_mass * eta**-0.6
    root = np.sqrt(max(1.0 - 4.0 * eta, 0.0))
    return total * (1 + root) / 2, total * (1 - root) / 2


def inspiral_duration(chirp_mass: float, f_low: float) -> float:
    """Leading-order time (s) a chirp spends above f_low before coalescence."""
    mc_s = chirp_mass * MSUN_S
    return 5.0 / 256.0 * mc_s ** (-5.0 / 3.0) * (np.pi * f_low) ** (-8.0 / 3.0)


def default_bank() -> list[tuple[float, float]]:
    """Coarse bank around stellar-mass binary black holes (m1, m2 in Msun)."""
    return [
        chirp_to_masses(float(mc), eta)
        for mc in np.arange(15.0, 45.1, 2.5)
        for eta in (0.18, 0.21, 0.23, 0.25)
    ]


def taylorf2(
    freqs: np.ndarray, m1: float, m2: float, f_low: float = 20.0, f_cut_factor: float = 1.0
) -> np.ndarray:
    """Inspiral template h~(f) ~ f^(-7/6) exp(-i psi(f)) with 2PN phase (arbitrary amplitude).

    Uses NumPy's FFT sign convention (exp(-2 pi i f t)). The template is zero outside
    [f_low, f_cut_factor * f_isco]. Values of f_cut_factor > 1 extrapolate the PN expansion
    beyond its validity and are meant for experimentation only.
    """
    mtot = (m1 + m2) * MSUN_S
    eta = m1 * m2 / (m1 + m2) ** 2
    f_cut = f_cut_factor / (6**1.5 * np.pi * mtot)
    band = (freqs >= f_low) & (freqs <= f_cut)
    h = np.zeros(len(freqs), dtype=complex)
    if not band.any():
        return h
    f = freqs[band]
    v = (np.pi * mtot * f) ** (1 / 3)
    psi = (3.0 / (128.0 * eta * v**5)) * (
        1.0
        + (20.0 / 9.0) * (743.0 / 336.0 + 11.0 * eta / 4.0) * v**2
        - 16.0 * np.pi * v**3
        + 10.0 * (3058673.0 / 1016064.0 + 5429.0 * eta / 1008.0 + 617.0 * eta**2 / 144.0) * v**4
    )
    h[band] = f ** (-7.0 / 6.0) * np.exp(-1j * (psi - np.pi / 4))
    return h


@dataclass
class MatchedFilterResult:
    snr: float
    peak_time_offset_s: float
    m1: float
    m2: float
    chirp_mass: float
    background_max: float
    background_median: float
    false_alarm_lower_bound: float
    n_templates: int
    n_off_windows: int
    times: np.ndarray = field(repr=False)
    rho: np.ndarray = field(repr=False)
    event_gps: float = 0.0

    def summary(self) -> dict:
        return {
            "snr": float(self.snr),
            "peak_time_offset_s": float(self.peak_time_offset_s),
            "m1_msun": float(self.m1),
            "m2_msun": float(self.m2),
            "chirp_mass_msun": float(self.chirp_mass),
            "background_max_snr": float(self.background_max),
            "background_median_snr": float(self.background_median),
            "false_alarm_probability_lower_bound": float(self.false_alarm_lower_bound),
            "n_templates": int(self.n_templates),
            "n_off_source_windows": int(self.n_off_windows),
        }


def matched_filter_search(
    times: np.ndarray,
    strain: np.ndarray,
    event_gps: float,
    bank: list[tuple[float, float]] | None = None,
    f_low: float = 20.0,
    window_s: float = 0.4,
    exclude_s: float = 1.0,  # seconds kept out of the background *before* the event
    f_cut_factor: float = 1.0,
) -> MatchedFilterResult:
    """Phase-maximised matched-filter SNR around `event_gps`, with an off-source background.

    SNR(t) = |z(t)| / sigma, z(t) = 4 * integral d~(f) h~*(f) / S(f) exp(2 pi i f t) df,
    sigma^2 = 4 * integral |h~(f)|^2 / S(f) df, S(f) = one-sided PSD (median Welch).
    The statistic is maximised over the bank; the on-source value is the peak within
    +-window_s/2 of the event. The background is the same peak taken over other windows
    of equal width in the same data (no signal), giving an empirical false-alarm bound.
    """
    times = np.asarray(times, dtype=float)
    n = len(strain)
    if n < 2:
        raise ValueError("Insufficient data points.")
    dt = float(times[1] - times[0])
    if n * dt > MAX_DURATION_S:
        raise ValueError(f"Data span {n * dt:.0f} s too long for matched filtering; crop it first.")
    fs = 1.0 / dt
    bank = bank or default_bank()

    x = (strain - np.mean(strain)) * tukey(n, alpha=min(1.0, 2.0 * fs / n))
    wf, wpsd = welch(x, fs=fs, nperseg=min(int(4 * fs), n), average="median")
    freqs = np.fft.rfftfreq(n, d=dt)
    psd = np.interp(freqs, wf, wpsd)
    data_f = dt * np.fft.rfft(x)
    df = 1.0 / (n * dt)

    rho_max = np.zeros(n)
    best_idx = np.zeros(n, dtype=int)
    for k, (m1, m2) in enumerate(bank):
        h = taylorf2(freqs, m1, m2, f_low, f_cut_factor)
        band = np.flatnonzero(h)
        if band.size == 0:
            continue
        sigma = np.sqrt(4.0 * df * np.sum(np.abs(h[band]) ** 2 / psd[band]))
        spec = np.zeros(n, dtype=complex)
        spec[band] = data_f[band] * np.conj(h[band]) / psd[band]
        rho = np.abs(4.0 * df * n * np.fft.ifft(spec)) / sigma
        better = rho > rho_max
        rho_max[better] = rho[better]
        best_idx[better] = k

    rel = times - event_gps
    on = np.abs(rel) <= window_s / 2
    if not on.any():
        raise ValueError("Event time falls outside the data span.")
    j_on = np.flatnonzero(on)[np.argmax(rho_max[on])]

    n_win = max(1, round(window_s * fs))
    edge = round(2 * fs)
    # A template output at time t depends on the data up to the longest template earlier, so a
    # loud signal leaks into later windows: keep them out of the background.
    tau_max = max(
        inspiral_duration((m1 * m2) ** 0.6 / (m1 + m2) ** 0.2, f_low) for m1, m2 in bank
    )
    skip_before, skip_after = max(exclude_s, window_s), tau_max + window_s + 0.5
    bg = []
    for s in range(edge, n - n_win - edge, n_win):
        start, end = rel[s], rel[s + n_win - 1]
        if end > -skip_before and start < skip_after:
            continue
        bg.append(float(rho_max[s : s + n_win].max()))
    if len(bg) < 5:
        raise ValueError("Not enough off-source data for a background estimate.")
    bg_arr = np.asarray(bg)
    snr = float(rho_max[j_on])
    m1, m2 = bank[best_idx[j_on]]
    mc = (m1 * m2) ** 0.6 / (m1 + m2) ** 0.2
    show = np.abs(rel) <= 1.0
    return MatchedFilterResult(
        snr=snr,
        peak_time_offset_s=float(rel[j_on]),
        m1=float(m1),
        m2=float(m2),
        chirp_mass=float(mc),
        background_max=float(bg_arr.max()),
        background_median=float(np.median(bg_arr)),
        false_alarm_lower_bound=float((1 + np.sum(bg_arr >= snr)) / (1 + len(bg_arr))),
        n_templates=len(bank),
        n_off_windows=len(bg_arr),
        times=rel[show],
        rho=rho_max[show],
        event_gps=float(event_gps),
    )
