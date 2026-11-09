# Simulación Monte Carlo de cadenas de desintegración radiactiva

Simulación estocástica de cadenas de decaimiento radiactivo en R usando el método de Monte Carlo.

**Autores:** Matías Berrío y Samuel Vargas

## Contenido

| Archivo | Descripción |
|---|---|
| `MC_validacion.R` | Valida el método comparando la simulación Monte Carlo con la solución analítica de Bateman para la cadena **Bi-212 → Tl-208 → Pb-208**. |
| `MC_aplicacion.R` | Aplica el método a la cadena del radón: **Rn-222 → Po-218 → Pb-214 → Bi-214 → Po-214 → Pb-210**. |

## Método

En cada paso de tiempo `dt`, el número de núcleos que decaen de cada isótopo se muestrea de una distribución binomial:

```
d ~ Binomial(N, p),   p = 1 - exp(-λ·dt),   λ = ln(2) / T½
```

Los núcleos que decaen pasan al siguiente isótopo de la cadena. Se repite la simulación varias veces y se promedian las trayectorias punto a punto.

## Requisitos

- R ≥ 4.1 (usa el operador `|>`)
- Paquetes: `ggplot2`, `tidyr`, `patchwork`

```r
install.packages(c("ggplot2", "tidyr", "patchwork"))
```

## Uso

Abre el proyecto en RStudio (`montecarlo-decaimiento.Rproj`) y ejecuta:

```r
source("MC_validacion.R")   # validación contra Bateman
source("MC_aplicacion.R")   # cadena del Rn-222
```

## Parámetros principales

| Parámetro | Validación | Aplicación |
|---|---|---|
| N₀ (átomos iniciales) | 1000 | 1000 |
| Tiempo máximo | 600 min | 200 h |
| Paso `dt` | 0.5 min | 0.5 h |
| Corridas | 30 | 100 |
| Semilla | 42 | 42 |
