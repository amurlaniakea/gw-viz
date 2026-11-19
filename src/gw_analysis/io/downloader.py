"""Download and load GWOSC strain files (HDF5)."""
from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
import requests
from gwosc import datasets
from gwosc.locate import get_event_urls

DETECTORS = ("H1", "L1", "V1")

# Fallback if the GWOSC catalogue lookup is unavailable (GPS of the merger).
KNOWN_EVENTS_GPS = {
    "GW150914": 1126259462.4,
    "GW170817": 1187008882.4,
    "GW190521": 1242442967.4,
}


@dataclass(frozen=True)
class Strain:
    """Strain time series with absolute GPS timestamps."""

    times: np.ndarray  # GPS seconds
    strain: np.ndarray
    sample_rate: float


def crop_strain(s: Strain, center_gps: float, span_s: float) -> Strain:
    """Keep `span_s` seconds centred on `center_gps` (whole series if span_s <= 0)."""
    if span_s <= 0:
        return s
    keep = np.abs(s.times - center_gps) <= span_s / 2
    if keep.sum() < 2:
        raise ValueError("Event time falls outside the data span.")
    return Strain(times=s.times[keep], strain=s.strain[keep], sample_rate=s.sample_rate)


def data_dir() -> Path:
    return Path(os.getenv("GW_DATA_DIR", "./data/raw"))


def _retrying_session() -> requests.Session:
    """HTTP session that retries dropped connections and 429/5xx answers with backoff."""
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry

    session = requests.Session()
    retry = Retry(total=4, backoff_factor=1.0, status_forcelist=(429, 500, 502, 503, 504),
                  allowed_methods=frozenset({"GET"}))
    session.mount("https://", HTTPAdapter(max_retries=retry))
    return session


def list_events(include_marginal: bool = False) -> list[str]:
    """Event names in the GWOSC catalogue (needs network); version suffixes are dropped.

    By default only names starting with ``GW`` are returned; the catalogue also holds marginal
    candidates named by date only (e.g. ``151008``), which are not the usual teaching examples.
    """
    names = datasets.find_datasets(type="event", session=_retrying_session())
    clean = {re.sub(r"-v\d+$", "", n) for n in names}
    if not include_marginal:
        clean = {n for n in clean if re.match(r"GW\d{6}", n)}
    return sorted(clean)


def event_gps(event_id: str) -> float:
    """GPS time of an event: GWOSC catalogue first, built-in table as fallback."""
    try:
        return float(datasets.event_gps(event_id))
    except Exception:  # noqa: BLE001 - any lookup failure falls back to the table
        if event_id in KNOWN_EVENTS_GPS:
            return KNOWN_EVENTS_GPS[event_id]
        raise ValueError(f"Unknown event {event_id!r} and GWOSC lookup failed.") from None


def _pick_url(urls: list[str]) -> str:
    # Prefer the short event file (32 s) over the 4096 s one when both exist.
    short = [u for u in urls if "-32." in u]
    return (short or urls)[0]


def download_event_strain(
    event_id: str, detector: str = "H1", dest_dir: Path | None = None
) -> Path:
    """Download the strain file for an event, with an atomic local cache."""
    if detector not in DETECTORS:
        raise ValueError(f"Detector must be one of {DETECTORS}, got {detector!r}.")
    dest_dir = Path(dest_dir) if dest_dir else data_dir()
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / f"{event_id}_{detector}.hdf5"
    if dest.exists():
        return dest

    urls = get_event_urls(event_id, detector=detector, sample_rate=4096, format="hdf5")
    if not urls:
        raise RuntimeError(f"GWOSC returned no 4096 Hz HDF5 files for {event_id}/{detector}.")
    url = _pick_url(list(urls))

    tmp = dest.with_suffix(".hdf5.part")  # never leave a truncated file looking like a cache hit
    try:
        with requests.get(url, stream=True, timeout=60) as r:
            r.raise_for_status()
            with open(tmp, "wb") as f:
                f.writelines(r.iter_content(1 << 16))
        os.replace(tmp, dest)
    finally:
        tmp.unlink(missing_ok=True)
    return dest


def _find_gwpy_dataset(f: h5py.File) -> h5py.Dataset:
    found: list[h5py.Dataset] = []

    def visit(_name: str, obj: object) -> None:
        if isinstance(obj, h5py.Dataset) and "x0" in obj.attrs and "dx" in obj.attrs:
            found.append(obj)

    f.visititems(visit)
    if not found:
        raise ValueError("Unrecognised HDF5 layout: no GWOSC 'strain/Strain' and no gwpy dataset.")
    return found[0]


def load_hdf5_strain(path: Path) -> Strain:
    """Read strain from a GWOSC-native HDF5 (strain/Strain + meta/GPSstart) or from a
    file written by gwpy's TimeSeries.write (single dataset with x0/dx attributes)."""
    with h5py.File(path, "r") as f:
        if "strain/Strain" in f:
            ds = f["strain/Strain"]
            dt = float(ds.attrs["Xspacing"])
            t0 = float(f["meta/GPSstart"][()])
        else:
            ds = _find_gwpy_dataset(f)
            dt = float(ds.attrs["dx"])
            t0 = float(ds.attrs["x0"])
        strain = ds[:].astype(float)
    times = t0 + np.arange(len(strain)) * dt
    return Strain(times=times, strain=strain, sample_rate=1.0 / dt)
