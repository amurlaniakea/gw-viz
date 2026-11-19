"""Programmatic equivalent of `gw-analyze GW150914`. Needs network the first time."""
import matplotlib

matplotlib.use("Agg")

from gw_analysis.io.downloader import download_event_strain, event_gps, load_hdf5_strain
from gw_analysis.processing.snr_calculator import analyze_signal
from gw_analysis.viz.plots import plot_full_analysis

path = download_event_strain("GW150914", "H1")
s = load_hdf5_strain(path)
m = analyze_signal(s.times, s.strain, event_gps=event_gps("GW150914"))
print(m)
print(plot_full_analysis(s.times, s.strain, m, "GW150914"))
