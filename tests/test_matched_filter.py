import numpy as np
import pytest
from scipy.signal import hilbert

from gw_analysis.processing.matched_filter import (
    chirp_to_masses,
    default_bank,
    matched_filter_search,
    taylorf2,
)
from tests.conftest import EVENT, FS, GPS0, inject_taylorf2, make_strain


def _times(n):
    return GPS0 + np.arange(n) / FS


def test_template_frequency_increases_with_time():
    n = int(32 * FS)
    f = np.fft.rfftfreq(n, 1 / FS)
    m1, m2 = chirp_to_masses(25.0, 0.23)
    h = taylorf2(f, m1, m2)
    sig = np.fft.irfft(h * np.exp(-2j * np.pi * f * 16.0), n=n)
    phase = np.unwrap(np.angle(hilbert(sig)))
    inst = np.gradient(phase) * FS / (2 * np.pi)
    t = np.arange(n) / FS - 16.0
    early = inst[np.argmin(np.abs(t + 0.30))]
    late = inst[np.argmin(np.abs(t + 0.05))]
    assert late > early  # a chirp: frequency rises toward coalescence


def test_masses_roundtrip():
    m1, m2 = chirp_to_masses(30.0, 0.23)
    assert (m1 * m2) ** 0.6 / (m1 + m2) ** 0.2 == pytest.approx(30.0)
    assert m1 > m2 > 0


def test_injection_is_recovered():
    noise = make_strain(False, seed=11)
    m1, m2 = chirp_to_masses(30.0, 0.23)
    data = inject_taylorf2(noise, m1, m2, rho_opt=25.0, t_c=16.0)
    r = matched_filter_search(_times(len(data)), data, EVENT)
    assert 18.0 < r.snr < 30.0
    assert abs(r.chirp_mass - 30.0) <= 2.6
    assert abs(r.peak_time_offset_s) < 0.05
    assert r.false_alarm_lower_bound == pytest.approx(1.0 / (1 + r.n_off_windows))


def test_noise_only_has_no_loud_peak():
    for seed in (1, 2):
        noise = make_strain(False, seed=seed)
        r = matched_filter_search(_times(len(noise)), noise, EVENT)
        assert r.snr < 7.5
        assert r.false_alarm_lower_bound > 0.01


def test_too_long_data_is_rejected():
    n = int(400 * FS)
    with pytest.raises(ValueError):
        matched_filter_search(_times(n), np.zeros(n), GPS0 + 200)


def test_bank_is_sane():
    bank = default_bank()
    assert len(bank) >= 40 and all(m1 >= m2 > 0 for m1, m2 in bank)


def test_loud_signal_does_not_pollute_background():
    noise = make_strain(False, seed=7, dur=48)
    m1, m2 = chirp_to_masses(25.0, 0.21)
    data = inject_taylorf2(noise, m1, m2, rho_opt=40.0, t_c=24.0)
    r = matched_filter_search(_times(len(data)), data, GPS0 + 24.0)
    assert r.snr > 30 and r.background_max < 8
