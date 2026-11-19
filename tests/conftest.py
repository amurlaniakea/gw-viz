import h5py
import numpy as np
import pytest

FS = 4096.0
GPS0 = 1126259446.0  # file start; event 16 s in, like GWOSC 32 s event files
EVENT = GPS0 + 16.0


def make_strain(inject: bool, seed: int = 0, dur: int = 32) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n = int(dur * FS)
    # coloured noise: white noise with a 1/f-ish low-frequency rise
    spec = np.fft.rfft(rng.normal(size=n))
    f = np.fft.rfftfreq(n, 1 / FS)
    spec /= np.sqrt(1 + (30.0 / np.maximum(f, 1.0)) ** 4)
    x = np.fft.irfft(spec, n=n) * 1e-21
    if inject:
        t = np.arange(n) / FS - 16.0
        m = (t > -0.2) & (t < 0.02)
        tau = 0.02 - t[m]
        freq = 35 * (1 + 8 * (0.22 / (tau + 0.02)) ** 0.375 / 8)
        x[m] += 8e-21 * np.sin(2 * np.pi * np.cumsum(freq) / FS) * np.hanning(m.sum())
    return x


@pytest.fixture
def noise():
    return make_strain(False)


@pytest.fixture
def signal():
    return make_strain(True)


@pytest.fixture
def gwosc_file(tmp_path, signal):
    p = tmp_path / "GW150914_H1.hdf5"
    with h5py.File(p, "w") as f:
        ds = f.create_dataset("strain/Strain", data=signal)
        ds.attrs["Xspacing"] = 1 / FS
        f.create_dataset("meta/GPSstart", data=int(GPS0))
    return p


def colored_psd(freqs):
    """Analytic one-sided PSD of make_strain() noise."""
    return (2.0 / FS) * 1e-42 / (1 + (30.0 / np.maximum(freqs, 1.0)) ** 4)


def inject_taylorf2(noise, m1, m2, rho_opt, t_c, f_low=20.0):
    """Add a TaylorF2 inspiral with known optimal SNR (given analytic PSD) ending at t_c seconds."""
    from gw_analysis.processing.matched_filter import taylorf2

    n = len(noise)
    f = np.fft.rfftfreq(n, 1 / FS)
    h = taylorf2(f, m1, m2, f_low)
    band = np.flatnonzero(h)
    sigma_unit = np.sqrt(4 * (1 / (n / FS)) * np.sum(np.abs(h[band]) ** 2 / colored_psd(f)[band]))
    amp = rho_opt / sigma_unit
    return noise + np.fft.irfft(amp * h * np.exp(-2j * np.pi * f * t_c) * FS, n=n)


@pytest.fixture(autouse=True)
def no_network(monkeypatch):
    """Tests must be hermetic: any real HTTP request fails immediately."""

    def blocked(*_a, **_k):
        raise RuntimeError("network access is blocked in tests")

    monkeypatch.setattr("requests.sessions.Session.request", blocked)
