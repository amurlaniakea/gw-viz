# Guía de prácticas (nivel grado / máster)

**Requisitos**: Python básico, nociones de Fourier y estadística. **Duración**: 4 sesiones de ~1 h (se pueden agrupar).
**Preparación del profesor**: instalar el paquete (README) y descargar GW150914 una vez; compartir `data/raw/` evita que 30 alumnos descarguen a la vez.

## Práctica 1 — La señal cruda no se ve
Ejecuta `gw-analyze GW150914 --detector H1` y abre `output/GW150914_report.png`.
- ¿Por qué el strain crudo parece ruido? ¿Qué orden de magnitud tiene (≈10⁻²¹)?
- Compara el panel filtrado con el blanqueado. ¿Qué aporta cada paso?
*Objetivo*: entender ruido coloreado, filtrado y blanqueado.

## Práctica 2 — Ruido y PSD
En el notebook, representa la ASD del panel 3. Localiza la banda sensible y las líneas estrechas (60 Hz y armónicos).
- ¿Por qué el ruido sube a baja frecuencia? ¿Y a alta frecuencia?
- ¿Qué pasaría con la amplitud filtrada si se aplicaran filtros notch?
*Objetivo*: leer una densidad espectral; relación señal-ruido en frecuencia.

## Práctica 3 — ¿Es robusta una estadística?
Ejecuta con `--span 32`, `64` y `128` y tabula `excess_power_sigma` y `peak_whitened_sigma`.
- ¿Cuánto varían? ¿Por qué usamos una mediana en la PSD (compara con `average="mean"` en `whiten_bandpass`)?
- Repite con L1. ¿Qué diferencia de tiempo de llegada esperas entre H1 y L1?
*Objetivo*: sensibilidad al método, estimadores robustos, reproducibilidad.

## Práctica 4 — Filtro adaptado y falsa alarma
Ejecuta con `--matched-filter` y abre `output/GW150914_snr.png`.
- ¿Qué masa chirp recupera? ¿Por qué el SNR es menor que el publicado (~24)? (pista: plantillas cortadas en el ISCO)
- Interpreta la cota de falsa alarma: ¿qué significa «(1 + k)/(1 + N)»? ¿Cómo la mejorarías con más ventanas?
- Voluntario: inyecta una señal con `tests/conftest.py::inject_taylorf2` en ruido sintético y comprueba que recuperas el SNR inyectado.
*Objetivo*: idea de filtro óptimo, fondo empírico y significancia estadística.

## Discusión final
1. ¿Qué diferencia hay entre «ver» una señal y «detectarla»?
2. ¿Por qué se exige coincidencia entre detectores en el análisis real?
3. ¿Qué añadiría un banco con fusión y ringdown (p. ej. IMRPhenom)?

## Cuidado con los eventos
Esta herramienta está calibrada para agujeros negros estelares (banda 35–350 Hz). Para sistemas de estrellas de neutrones, como GW170817,
la señal dura minutos y se concentra en otras frecuencias: sirve como ejemplo de **límite del método**, no como práctica de detección.
