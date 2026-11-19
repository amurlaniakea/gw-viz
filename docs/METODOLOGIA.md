# Metodología

Este documento describe qué hace exactamente el código, para que pueda auditarse y usarse en clase.

## 1. Datos

Archivos HDF5 de strain a 4096 Hz del [GWOSC](https://gwosc.org/). El cargador lee del propio archivo la frecuencia de
muestreo y el instante GPS inicial (formato nativo `strain/Strain` + `meta/GPSstart`, o el escrito por GWpy con atributos `x0`/`dx`).
El instante del evento sale del catálogo GWOSC (con tabla interna de respaldo). Los datos se recortan a `--span` segundos
(64 s por defecto) centrados en el evento, de modo que el resultado no depende de la longitud del archivo descargado.

## 2. Preprocesado

1. Se resta la media y se aplica una ventana de Tukey con taper fijo de 1 s en cada extremo.
2. PSD unilateral con el método de Welch (segmentos de 4 s, **promedio por mediana**). La mediana evita que la propia señal
   sesgue la estimación del ruido, que es lo que ocurría con el promedio por media.
3. Blanqueado en frecuencia (división por la ASD) y filtro Butterworth paso banda de orden 4, fase cero, 35–350 Hz.

## 3. Exceso de potencia

Ventana de 0,3 s (de −0,225 s a +0,075 s respecto al evento). Se calcula la energía de la señal blanqueada en esa ventana
y se compara con la de las demás ventanas del archivo, excluyendo 1 s alrededor del evento y 2 s en los bordes:

    excess_power_sigma = (E_on − media(E_off)) / desviación(E_off)

`peak_whitened_sigma` es el máximo de |señal blanqueada| en la ventana dividido por la desviación típica del ruido (de las ventanas fuera de la fuente).
`peak_amp_bandpassed` es el pico del strain filtrado 35–350 Hz, **sin** filtros notch: puede incluir líneas instrumentales (60 Hz y armónicos),
por lo que no es la amplitud astrofísica del evento.

## 4. Filtro adaptado (educativo)

Plantillas en el dominio de frecuencias, aproximación estacionaria de fase (TaylorF2) con fase hasta 2PN:

    h(f) ∝ f^(−7/6) · exp(−i ψ(f)),   ψ(f) = (3 / 128η) v⁻⁵ [1 + (20/9)(743/336 + 11η/4) v² − 16π v³ + 10(3058673/1016064 + 5429η/1008 + 617η²/144) v⁴]

con v = (π M f)^(1/3), cortada entre 20 Hz y la frecuencia del ISCO. El banco cubre masas chirp de 15 a 45 M☉ y η ∈ {0,18; 0,21; 0,23; 0,25}.

    z(t) = 4 ∫ d(f) h*(f) / S(f) · exp(2π i f t) df,   σ² = 4 ∫ |h(f)|² / S(f) df,   SNR(t) = |z(t)| / σ

El SNR se maximiza sobre la fase (módulo complejo) y sobre el banco. El valor «en fuente» es el máximo en ±0,2 s del evento.
El **fondo** son los máximos del mismo estadístico en ventanas de igual anchura sin señal; se excluyen las ventanas que la cola de las
plantillas más largas podría contaminar. La cota de falsa alarma es (1 + nº de ventanas ≥ SNR) / (1 + nº de ventanas).

## 5. Validación

- **Inyección y recuperación**: se inyecta una plantilla con SNR óptimo conocido (calculado con la PSD analítica del ruido sintético) en ruido
  gaussiano coloreado y se recupera. Con SNR inyectado de 10, 15, 25 y 40 se recuperó 8,7–11,2; 13,6–16,1; 23,5–25,8 y 38,1–40,5
  (3 realizaciones cada uno), con la masa chirp correcta. Solo ruido: máximo del banco 3,3–4,3.
- **Convención de signo**: un test comprueba que la frecuencia de la plantilla crece con el tiempo (un chirp) y no al revés.
- **Estabilidad**: con GW150914 (H1) la significancia de exceso de potencia fue 18,7σ, 18,1σ y 17,4σ para ventanas de 32, 64 y 128 s.
- **Datos reales (GW150914, `--span 64`)**: filtro adaptado con SNR 10,8 en H1 (fondo máximo 4,9) y 6,5 en L1 (fondo máximo 5,1); 138 ventanas fuera de la fuente en ambos casos. La suma en cuadratura (≈12,6) es ~la mitad del SNR de red publicado (~24), compatible con plantillas cortadas en el ISCO. Las masas del mejor ajuste (37 + 37 M☉) son indicativas, no una estimación de parámetros.
- **Cargadores y CLI**: tests con archivos HDF5 sintéticos en ambos formatos, descarga atómica y caché.

## 6. Limitaciones

Sin fusión ni ringdown (el SNR es una cota inferior), sin estimación de parámetros, sin coincidencia entre detectores, sin tratamiento de glitches,
sin filtros notch, banco grueso. Es una herramienta didáctica; para resultados científicos usa PyCBC, GWpy o LALSuite.

## 7. Referencias

- Abbott et al., «Observation of Gravitational Waves from a Binary Black Hole Merger», *Phys. Rev. Lett.* 116, 061102 (2016).
- Abbott et al., «Open data from the first and second observing runs of Advanced LIGO and Advanced Virgo», *SoftwareX* 13, 100658 (2021).
- Allen et al., «FINDCHIRP», *Phys. Rev. D* 85, 122006 (2012).
- Blanchet, «Gravitational Radiation from Post-Newtonian Sources and Inspiralling Compact Binaries», *Living Rev. Relativ.* 17, 2 (2014).
- Welch, *IEEE Trans. Audio Electroacoust.* 15, 70 (1967).
