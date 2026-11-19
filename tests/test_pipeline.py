import numpy as np
import pytest

from gw_analysis.io import downloader
from gw_analysis.io.downloader import load_hdf5_strain
from gw_analysis.processing.snr_calculator import analyze_signal
from tests.conftest import EVENT, FS, GPS0


def test_loader_reads_gps_and_rate(gwosc_file):
    s = load_hdf5_strain(gwosc_file)
    assert s.sample_rate == pytest.approx(FS)
    assert s.times[0] == pytest.approx(GPS0)
    assert s.times[-1] - s.times[0] == pytest.approx(32 - 1 / FS)


def test_injected_signal_is_significant_and_noise_is_not(gwosc_file, noise):
    s = load_hdf5_strain(gwosc_file)
    on = analyze_signal(s.times, s.strain, event_gps=EVENT)
    off = analyze_signal(s.times, noise, event_gps=EVENT)
    assert on["excess_power_sigma"] > 8
    assert abs(off["excess_power_sigma"]) < 4
    assert 30 < on["central_freq_hz"] < 350
    assert on["peak_amp_bandpassed"] < np.max(np.abs(s.strain)) * 10


def test_event_outside_data_raises(gwosc_file):
    s = load_hdf5_strain(gwosc_file)
    with pytest.raises(ValueError):
        analyze_signal(s.times, s.strain, event_gps=GPS0 + 1000)


def test_download_is_atomic_and_cached(tmp_path, monkeypatch):
    monkeypatch.setattr(downloader, "get_event_urls", lambda *a, **k: ["http://x/a-4096.hdf5", "http://x/a-32.hdf5"])
    calls = []

    class R:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def raise_for_status(self): pass
        def iter_content(self, n): yield b"abc"

    monkeypatch.setattr(downloader.requests, "get", lambda url, **k: calls.append(url) or R())
    p = downloader.download_event_strain("GW150914", "H1", tmp_path)
    assert p.read_bytes() == b"abc" and calls == ["http://x/a-32.hdf5"]
    assert not list(tmp_path.glob("*.part"))
    downloader.download_event_strain("GW150914", "H1", tmp_path)
    assert len(calls) == 1  # cache hit


def test_failed_download_leaves_no_cache(tmp_path, monkeypatch):
    monkeypatch.setattr(downloader, "get_event_urls", lambda *a, **k: ["http://x/a.hdf5"])

    def boom(*a, **k):
        raise OSError("net down")

    monkeypatch.setattr(downloader.requests, "get", boom)
    with pytest.raises(OSError):
        downloader.download_event_strain("GW150914", "H1", tmp_path)
    assert not list(tmp_path.iterdir())


def test_bad_detector():
    with pytest.raises(ValueError):
        downloader.download_event_strain("GW150914", "X1")


def test_loader_reads_gwpy_written_file(tmp_path, signal):
    import h5py

    p = tmp_path / "gwpy.hdf5"
    with h5py.File(p, "w") as f:
        ds = f.create_dataset("H1:GWOSC-4KHZ_R1_STRAIN", data=signal)
        ds.attrs["dx"] = 1 / FS
        ds.attrs["x0"] = GPS0
    s = load_hdf5_strain(p)
    assert s.sample_rate == pytest.approx(FS) and s.times[0] == pytest.approx(GPS0)


def test_loader_rejects_unknown_layout(tmp_path):
    import h5py

    p = tmp_path / "bad.hdf5"
    with h5py.File(p, "w") as f:
        f.create_dataset("x", data=np.zeros(10))
    with pytest.raises(ValueError):
        load_hdf5_strain(p)


def test_crop_keeps_requested_span(gwosc_file):
    from gw_analysis.io.downloader import crop_strain

    s = load_hdf5_strain(gwosc_file)
    c = crop_strain(s, EVENT, 24)
    assert c.times[0] >= EVENT - 12 and c.times[-1] <= EVENT + 12
    a = analyze_signal(c.times, c.strain, event_gps=EVENT)
    assert a["excess_power_sigma"] > 8
    with pytest.raises(ValueError):
        crop_strain(s, GPS0 + 5000, 24)
    assert crop_strain(s, EVENT, 0) is s


def test_significance_is_stable_across_analysed_span(signal):
    """Median-averaged PSD keeps the injected signal from biasing its own whitening."""
    from gw_analysis.io.downloader import Strain, crop_strain

    t = GPS0 + np.arange(len(signal)) / FS
    s = Strain(t, signal, FS)
    vals = [
        analyze_signal(c.times, c.strain, event_gps=EVENT)["excess_power_sigma"]
        for c in (crop_strain(s, EVENT, span) for span in (0, 24))
    ]
    assert max(vals) / min(vals) < 1.25


def test_list_events_strips_versions(monkeypatch):
    monkeypatch.setattr(downloader.datasets, "find_datasets", lambda **k: ["GW150914-v3", "GW150914-v1", "GW170817-v3"])
    assert downloader.list_events() == ["GW150914", "GW170817"]


def test_cli_json_and_errors(tmp_path, gwosc_file, monkeypatch, capsys):
    import json

    from gw_analysis import cli

    monkeypatch.setattr(cli, "event_gps", lambda _e: EVENT)
    code = cli.main(["GW150914", "--data-dir", str(gwosc_file.parent), "--output-dir", str(tmp_path),
                     "--json", "--no-plot", "--span", "24"])
    out = json.loads(capsys.readouterr().out)
    assert code == 0 and out["event"] == "GW150914" and out["excess_power_sigma"] > 8
    # no cached file and no network -> clean error, exit code 1 (no real download attempted)
    assert cli.main(["GW150914", "--detector", "L1", "--data-dir", str(tmp_path / "none")]) == 1
    assert "Error" in capsys.readouterr().err
    with pytest.raises(SystemExit):
        cli.main([])


def test_list_events_falls_back_when_catalogue_is_down(monkeypatch, capsys):
    from gw_analysis import cli

    def down():
        raise ConnectionError("catalogue down")

    monkeypatch.setattr(cli, "list_events", down)
    assert cli.main(["--list-events"]) == 0
    cap = capsys.readouterr()
    assert "GW150914" in cap.out and "Warning" in cap.err


def test_retrying_session_retries_get():
    s = downloader._retrying_session()
    retry = s.get_adapter("https://gwosc.org").max_retries
    assert retry.total == 4 and 503 in retry.status_forcelist
