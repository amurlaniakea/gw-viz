"""Command line entry point: `gw-analyze GW150914 --detector H1`."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from gw_analysis.io.downloader import (
    crop_strain,
    download_event_strain,
    event_gps,
    load_hdf5_strain,
)
from gw_analysis.processing.snr_calculator import analyze_signal
from gw_analysis.viz.plots import plot_full_analysis


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="gw-analyze", description="Analyze a GWOSC event.")
    p.add_argument("event_id", help="Event name, e.g. GW150914")
    p.add_argument("--detector", default="H1", choices=["H1", "L1", "V1"])
    p.add_argument("--data-dir", type=Path, default=None, help="Cache dir (default ./data/raw)")
    p.add_argument("--span", type=float, default=64.0,
                   help="Seconds of data centred on the event to analyse (0 = whole file)")
    p.add_argument("--output-dir", type=Path, default=Path("output"))
    a = p.parse_args(argv)

    try:
        path = download_event_strain(a.event_id, a.detector, a.data_dir)
        s = load_hdf5_strain(path)
        gps = event_gps(a.event_id)
        s = crop_strain(s, gps, a.span)  # same result whether the file is 32 s, 60 s or 4096 s
        m = analyze_signal(s.times, s.strain, event_gps=gps)
        png = plot_full_analysis(s.times, s.strain, m, a.event_id, a.output_dir)
    except Exception as e:  # noqa: BLE001 - report cleanly, non-zero exit
        print(f"Error: {e}", file=sys.stderr)
        return 1

    print(f"{a.event_id} ({a.detector}): {len(s.strain)} samples @ {m['sample_rate_hz']:.0f} Hz")
    print(f"  excess-power significance : {m['excess_power_sigma']:.1f} sigma (NOT matched-filter SNR)")
    print(f"  peak whitened amplitude   : {m['peak_whitened_sigma']:.1f} sigma")
    print(f"  peak band-passed strain   : {m['peak_amp_bandpassed']:.2e}")
    print(f"  central frequency         : {m['central_freq_hz']:.0f} Hz")
    print(f"  report                    : {png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
