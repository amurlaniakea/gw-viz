"""Download and load GWOSC strain files (HDF5)."""
from __future__ import annotations

import os
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


def data_dir() -> Path:
    return Path(os.getenv("GW_DATA_DIR", "./data/raw"))


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
