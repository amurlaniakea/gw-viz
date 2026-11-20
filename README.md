# GW Analysis Toolkit

[![ci](https://github.com/amurlaniakea/gw-viz/actions/workflows/ci.yml/badge.svg)](https://github.com/amurlaniakea/gw-viz/actions/workflows/ci.yml)
![licencia](https://img.shields.io/badge/licencia-AGPL--3.0--or--later-blue)

Herramienta en Python, legible y reproducible, para analizar datos abiertos de ondas gravitacionales
(LIGO/Virgo/KAGRA) del [GWOSC](https://gwosc.org/). Pensada para **consulta rápida** (línea de comandos y salida JSON)
y para **docencia universitaria** (código corto, tests, notebook y guía de prácticas).

> *«Tecnosfera al servicio de la biosfera.»*

*English summary:* a small, readable and tested Python toolkit that downloads open LIGO/Virgo strain data,
whitens and band-passes it, measures excess-power significance against off-source data, and runs an
educational matched filter (2PN TaylorF2 inspiral bank) validated by injection/recovery. Not a replacement for
LALSuite/PyCBC; see *Limitaciones*.

**Demo y resultados en la web:** <https://amurlaniakea.github.io/gw-viz/> — instantáneas reales de GW150914 (H1/L1) y una demo interactiva
con **datos sintéticos** que muestra cómo un filtro adaptado recupera una señal inyectada. La herramienta en sí se ejecuta en Python, no en el navegador.

## Instalación

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

## Uso

```bash
gw-analyze GW150914 --detector H1                      # análisis básico + figura
gw-analyze GW150914 --detector H1 --matched-filter     # + matched filter educativo
gw-analyze GW150914 --span 32 --json --no-plot         # salida máquina-legible
gw-analyze --list-events                               # catálogo GWOSC, eventos GW… (requiere red)
gw-analyze --list-events --all-events                  # incluye candidatos marginales (nombres sin «GW»)
```

Salida real con GW150914 (H1, `--span 64`):

```
GW150914 (H1): 262144 samples @ 4096 Hz
  excess-power significance : 18.1 sigma (NOT matched-filter SNR)
  peak whitened amplitude   : 8.1 sigma
  peak band-passed strain   : 2.43e-21
  central frequency         : 57 Hz
```

Con `--matched-filter` (mismos datos, `--span 64`):

| Detector | SNR filtro adaptado | Mayor SNR en el fondo | Mejor plantilla (indicativa) | Cota de falsa alarma |
|---|---|---|---|---|
| H1 | 10.8 | 4.9 | m₁ = m₂ = 37 M☉ (M_chirp ≈ 32.5 M☉) | ≤ 0.007 (138 ventanas) |
| L1 | 6.5 | 5.1 | m₁ = m₂ = 37 M☉ (M_chirp ≈ 32.5 M☉) | ≤ 0.007 (138 ventanas) |

La suma en cuadratura (√(10.8² + 6.5²) ≈ 12.6) es aproximadamente la mitad del SNR de red publicado (~24), como cabe esperar
de plantillas solo de inspiral. En L1 la separación respecto al fondo (6.5 frente a 5.1) es modesta.

La primera ejecución descarga el archivo a `data/raw/` (o a `$GW_DATA_DIR`); las siguientes usan la caché.
Con `--span 32`, `64` y `128` la significancia sale 18.7σ, 18.1σ y 17.4σ: es estable frente a la ventana analizada.

## Qué mide (y qué no)

| Magnitud | Qué es |
|---|---|
| `excess_power_sigma` | Energía blanqueada de la ventana del evento frente a ventanas sin señal del mismo archivo (z-score). **No es SNR de filtro adaptado.** |
| `peak_whitened_sigma` | Pico de la señal blanqueada, en desviaciones típicas del ruido. |
| `peak_amp_bandpassed` | Pico del strain filtrado 35–350 Hz en la ventana (sin filtros notch). |
| `matched_filter.snr` | SNR de un filtro adaptado con plantillas de inspiral TaylorF2 (2PN) cortadas en el ISCO: **cota inferior**, no reproduce los SNR publicados (~24 de red en GW150914). |
| `false_alarm_probability_lower_bound` | Cota empírica: (1 + nº de ventanas sin señal con SNR ≥ el medido) / (1 + nº de ventanas). |

Detalles, fórmulas y validación en [`docs/METODOLOGIA.md`](docs/METODOLOGIA.md).

## Para clase

[`docs/AULA.md`](docs/AULA.md) propone cuatro prácticas (señal cruda vs. filtrada, ruido y PSD, robustez estadística,
filtro adaptado y probabilidad de falsa alarma) y preguntas de discusión. El notebook
[`notebooks/tutorial_gw150914.ipynb`](notebooks/tutorial_gw150914.ipynb) las recorre paso a paso.

## Estructura

- `src/gw_analysis/io/`: descarga con caché atómica y lectura de HDF5 (GPS y frecuencia de muestreo se leen del archivo).
- `src/gw_analysis/processing/`: blanqueado y métricas (`snr_calculator.py`), filtro adaptado (`matched_filter.py`).
- `src/gw_analysis/viz/`: figuras. `src/gw_analysis/cli.py`: comando `gw-analyze`.
- `docs/index.html` + `docs/demo.js`: página estática (GitHub Pages) con resultados reales y la demo interactiva; sin dependencias externas.
- `tests/`: 22 pruebas Python herméticas (la red está bloqueada) y 8 pruebas JS de la demo (paridad numérica con Python, inyección y recuperación).

```bash
pytest -q && ruff check . && node tests/demo.test.js
```

## Limitaciones

- Plantillas solo de inspiral (sin fusión ni ringdown): masas recuperadas **indicativas**, no estimación de parámetros.
- El banco cubre masas chirp de 15 a 45 M☉; sistemas más masivos saturan el borde del banco y el SNR queda aún más infravalorado.
- Un solo detector por ejecución; sin coincidencia H1/L1, sin tratamiento de glitches ni filtros notch.
- Validado con datos sintéticos y con GW150914 (H1). Otros eventos, en especial sistemas de estrellas de neutrones
  de larga duración como GW170817, requieren otra banda y otro banco de plantillas.
- Para análisis científicos usa [GWpy](https://gwpy.github.io/), [PyCBC](https://pycbc.org/) o LALSuite.

## Cómo citar

Consulta [`CITATION.cff`](CITATION.cff). Si usas datos de GWOSC, reconoce su origen:
*«This research has made use of data obtained from the Gravitational Wave Open Science Center (gwosc.org)»*.

Referencias: Abbott et al., *Phys. Rev. Lett.* 116, 061102 (2016) (GW150914); Abbott et al., *SoftwareX* 13, 100658 (2021) (datos abiertos O1/O2).

## Licencia

AGPL-3.0-or-later. © Pedro Sordo Martínez.
