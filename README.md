# Simulación Monte Carlo de cadenas de desintegración radiactiva

Simulación estocástica de cadenas de decaimiento radiactivo en R usando el método de Monte Carlo. El método es general: se valida con una cadena corta que tiene solución analítica y se aplica a la cadena del Rn-222 como caso de estudio, pero puede adaptarse a cualquier cadena lineal de decaimiento (ver [Escalabilidad](#escalabilidad)).

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

## Escalabilidad

Aunque `MC_aplicacion.R` está configurado para el Rn-222, nada del método depende de esa cadena en particular. Cada isótopo solo necesita su vida media, así que el mismo esquema sirve para otras cadenas naturales (U-238, Th-232, U-235), para cadenas artificiales o para cadenas de cualquier longitud.

### Cambiar de cadena

1. **Vidas medias:** reemplaza el vector `T12` con los isótopos de la nueva cadena, en orden y en una misma unidad de tiempo. Usa `Inf` para el isótopo estable final.
2. **Actualización de poblaciones:** las líneas que actualizan `N["Rn-222"]`, `N["Po-218"]`, … están escritas para esta cadena. Se pueden reemplazar por un ciclo general que funciona para cualquier cadena lineal:

   ```r
   # d[k] = núcleos del isótopo k que decayeron en este paso
   for (k in seq_len(n_iso)) {
     entrada <- if (k > 1) d[k - 1] else 0
     N[k] <- N[k] + entrada - d[k]
   }
   ```

3. **Gráficas:** ajusta el vector `colores` y los grupos de isótopos de cada panel (`intermedios`).
4. **Escala de tiempo:** elige `t_max` según la vida media más larga y `dt` según la más corta (ver abajo).

### Elección del paso de tiempo

El paso `dt` debe ser pequeño comparado con la vida media más corta que se quiera resolver. Si un isótopo tiene una vida media mucho menor que `dt` (como el Po-214, de 164 µs), puede tratarse como instantáneo, igual que se hace en este proyecto. Esto permite simular cadenas con escalas de tiempo muy distintas sin un costo computacional excesivo.

### Costo computacional

El tiempo de ejecución crece de forma aproximadamente lineal con:

- el número de isótopos de la cadena,
- el número de pasos (`t_max / dt`),
- el número de corridas (`n_runs`).

No depende de `N0`, porque se muestrea con `rbinom` el número total de decaimientos por paso en lugar de simular cada átomo. Por eso se pueden usar poblaciones iniciales grandes (por ejemplo, 10⁶ átomos) sin que la simulación se vuelva más lenta. Más corridas reducen el ruido del promedio, cuyo error disminuye como `1/√n_runs`.

### Posibles extensiones

- **Ramificaciones:** algunos isótopos decaen por dos caminos (por ejemplo, el Bi-214 produce Tl-210 en una fracción muy pequeña de los casos). Se puede incluir dividiendo los decaimientos con `rbinom` según la razón de ramificación.
- **Fuentes externas:** agregar un término de producción constante para modelar, por ejemplo, la entrada continua de radón a un recinto cerrado.
- **Intervalos de confianza:** guardar todas las corridas (como en `MC_validacion.R`) para graficar bandas de dispersión además del promedio.

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
