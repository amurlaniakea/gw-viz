"""Three-panel report: band-passed strain, whitened zoom on the event, ASD."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, sosfiltfilt, welch

from gw_analysis.processing.snr_calculator import BAND, whiten_bandpass


def plot_full_analysis(times, strain, metrics, event_id, output_dir="./output") -> Path:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    fs = metrics["sample_rate_hz"]
    t = np.asarray(times) - metrics["event_gps"]

    sos = butter(4, BAND, btype="bandpass", fs=fs, output="sos")
    bp = sosfiltfilt(sos, strain - np.mean(strain))
    w = whiten_bandpass(strain, fs)

    fig, ax = plt.subplots(3, 1, figsize=(10, 12))
    ax[0].plot(t, bp, color="#38bdf8", lw=0.8)
    ax[0].set(title=f"{event_id}: band-passed strain {BAND[0]:.0f}-{BAND[1]:.0f} Hz",
              xlabel="Time from event (s)", ylabel="Strain h")
    zoom = (t > -0.4) & (t < 0.2)
    ax[1].plot(t[zoom], w[zoom], color="#a78bfa", lw=1.2)
    ax[1].set(title="Whitened strain, zoom on the event",
              xlabel="Time from event (s)", ylabel="Whitened strain (arb.)")
    f, psd = welch(strain, fs=fs, nperseg=min(int(4 * fs), len(strain)))
    ax[2].loglog(f[1:], np.sqrt(psd[1:]), color="#22c55e", lw=1.2)
    for edge in BAND:
        ax[2].axvline(edge, color="red", ls=":")
    ax[2].set(title="Amplitude spectral density", xlabel="Frequency (Hz)",
              ylabel=r"ASD (strain/$\sqrt{\mathrm{Hz}}$)", xlim=(10, fs / 2))
    for a in ax:
        a.grid(True, alpha=0.3)
    fig.tight_layout()
    path = out / f"{event_id}_report.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path
