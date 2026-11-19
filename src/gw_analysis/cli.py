"""Command line entry point: `gw-analyze GW150914 --detector H1`."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from gw_analysis import __version__
from gw_analysis.io.downloader import (
    crop_strain,
    download_event_strain,
    event_gps,
    list_events,
    load_hdf5_strain,
)
from gw_analysis.processing.matched_filter import matched_filter_search
from gw_analysis.processing.snr_calculator import analyze_signal
from gw_analysis.viz.plots import plot_full_analysis, plot_snr_timeseries


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="gw-analyze", description="Analyze a GWOSC event.")
    p.add_argument("event_id", nargs="?", help="Event name, e.g. GW150914")
    p.add_argument("--detector", default="H1", choices=["H1", "L1", "V1"])
    p.add_argument("--data-dir", type=Path, default=None, help="Cache dir (default ./data/raw)")
    p.add_argument("--output-dir", type=Path, default=Path("output"))
    p.add_argument("--span", type=float, default=64.0,
                   help="Seconds of data centred on the event to analyse (0 = whole file)")
    p.add_argument("--matched-filter", action="store_true",
                   help="Also run the educational matched filter (inspiral-only bank)")
    p.add_argument("--no-plot", action="store_true", help="Skip figures")
    p.add_argument("--json", action="store_true", help="Print results as JSON only")
    p.add_argument("--list-events", action="store_true", help="List GWOSC events and exit")
    p.add_argument("--version", action="version", version=f"gw-analysis {__version__}")
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    a = parser.parse_args(argv)

    if a.list_events:
        try:
            print("\n".join(list_events()))
        except Exception as e:  # noqa: BLE001 - network/catalogue errors
            print(f"Error: {e}", file=sys.stderr)
            return 1
        return 0
    if not a.event_id:
        parser.error("event_id is required (or use --list-events)")

    try:
        path = download_event_strain(a.event_id, a.detector, a.data_dir)
        s = load_hdf5_strain(path)
        gps = event_gps(a.event_id)
        s = crop_strain(s, gps, a.span)  # same analysis window whatever the file length
        m = analyze_signal(s.times, s.strain, event_gps=gps)
        out: dict = {"event": a.event_id, "detector": a.detector, "samples": len(s.strain), **m}
        if not a.no_plot:
            out["report"] = str(plot_full_analysis(s.times, s.strain, m, a.event_id, a.output_dir))
        mf = None
        if a.matched_filter:
            mf = matched_filter_search(s.times, s.strain, gps)
            out["matched_filter"] = mf.summary()
            if not a.no_plot:
                out["snr_plot"] = str(plot_snr_timeseries(mf, a.event_id, a.output_dir))
    except Exception as e:  # noqa: BLE001 - report cleanly, non-zero exit
        print(f"Error: {e}", file=sys.stderr)
        return 1

    if a.json:
        print(json.dumps(out, indent=2))
        return 0
    print(f"{a.event_id} ({a.detector}): {out['samples']} samples @ {m['sample_rate_hz']:.0f} Hz")
    print(f"  excess-power significance : {m['excess_power_sigma']:.1f} sigma (NOT matched-filter SNR)")
    print(f"  peak whitened amplitude   : {m['peak_whitened_sigma']:.1f} sigma")
    print(f"  peak band-passed strain   : {m['peak_amp_bandpassed']:.2e}")
    print(f"  central frequency         : {m['central_freq_hz']:.0f} Hz")
    if mf:
        print("  --- matched filter (educational, inspiral-only templates: a lower bound) ---")
        print(f"  SNR                       : {mf.snr:.1f}  (best background window: {mf.background_max:.1f})")
        print(f"  best template             : m1={mf.m1:.0f} m2={mf.m2:.0f} Msun (Mc={mf.chirp_mass:.1f})")
        print(f"  false-alarm prob. (bound) : <= {mf.false_alarm_lower_bound:.3f} "
              f"({mf.n_off_windows} off-source windows)")
    for key in ("report", "snr_plot"):
        if key in out:
            print(f"  {key:<26}: {out[key]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
