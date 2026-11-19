# GW Analysis Toolkit

Análisis reproducible de datos abiertos de ondas gravitacionales (LIGO/Virgo/KAGRA) descargados de [GWOSC](https://www.gwosc.org/).

> *«Tecnosfera al servicio de la biosfera.»*

## Instalación

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

## Uso

```bash
gw-analyze GW150914 --detector H1
```

Qué hace:

1. Descarga el HDF5 del evento (4096 Hz) y lo cachea en `data/raw/` (o en `$GW_DATA_DIR`).
2. Blanquea y filtra en 35–350 Hz con una PSD de Welch estimada de los propios datos.
3. Mide la ventana del evento frente a ventanas fuera de la fuente del mismo archivo.
4. Escribe `output/GW150914_report.png` (señal filtrada, zoom blanqueado, ASD).

## Qué mide (y qué no)

Las métricas (`excess_power_sigma`, `peak_whitened_sigma`, `central_freq_hz`) son **significancia de exceso de potencia**
respecto al ruido del propio archivo. **No es SNR de filtro adaptado**: los SNR publicados
(p. ej. ~24 de red para GW150914) salen de matched filtering con bancos de plantillas y no coincidirán con estos números.

## Estructura

- `src/gw_analysis/io/`: descarga y lectura de HDF5 (tiempos GPS y frecuencia de muestreo leídos del propio archivo).
- `src/gw_analysis/processing/`: blanqueado, filtrado y métricas (NumPy/SciPy).
- `src/gw_analysis/viz/`: informe de 3 paneles.
- `src/gw_analysis/cli.py`: comando `gw-analyze`.
- `notebooks/`: tutorial paso a paso. `examples/`: script equivalente al CLI.
- `tests/`: pruebas con datos sintéticos (`pytest`), sin red.

## Desarrollo

```bash
pytest && ruff check .
```

## Licencia

AGPL-3.0-or-later. © Pedro Sordo Martínez.
