# Codigo de aplicacion del metodo de Monte Carlo

# Matias Berrio
# Samuel Vargas

library(ggplot2)
library(tidyr)

# Parámetros físicos
#  T½ en horas

T12 <- c(
  "Rn-222" = 3.8235 * 24,   # 3.8235 días → 91.76 h
  "Po-218" = 3.05 / 60,     # 3.05 min    → 0.05083 h
  "Pb-214" = 26.8 / 60,     # 26.8 min    → 0.4467 h
  "Bi-214" = 19.7 / 60,     # 19.7 min    → 0.3283 h
  "Po-214" = 0,             # instantáneo (164 µs)
  "Pb-210" = Inf            # estable
)

lambda <- log(2) / T12
lambda["Po-214"] <- Inf     # señal de decaimiento instantáneo
lambda["Pb-210"] <- 0       # estable

N0     <- 1000              # átomos iniciales de Rn-222
t_max  <- 200               # horas (~2 vidas medias del Rn-222)
dt     <- 1/2             # paso de tiempo [h] — fino para capturar hijos cortos
n_runs <- 100               # corridas Monte Carlo

tiempo  <- seq(0, t_max, by = dt)
n_steps <- length(tiempo)
nombres <- names(T12)
n_iso   <- length(nombres)

# Probabilidades de decaimiento por paso dt (Po-214 e Pb-210 tratados aparte)
p <- 1 - exp(-lambda * dt)
p["Po-214"] <- 1    # decae siempre (instantáneo)
p["Pb-210"] <- 0    # nunca decae

# Simulación Monte Carlo 

set.seed(42)

# Matriz de sumas acumuladas: fila = paso de tiempo, col = isótopo
suma <- matrix(0, nrow = n_steps, ncol = n_iso)
colnames(suma) <- nombres

for (r in seq_len(n_runs)) {
  
  # Estado inicial
  N <- setNames(rep(0, n_iso), nombres)
  N["Rn-222"] <- N0
  
  trayectoria <- matrix(0, nrow = n_steps, ncol = n_iso)
  colnames(trayectoria) <- nombres
  trayectoria[1, ] <- N
  
  for (i in 2:n_steps) {
    
    # Decaimientos de cada isótopo (excepto Pb-210 estable)
    d <- setNames(rep(0, n_iso), nombres)
    for (iso in nombres[-n_iso]) {          # excluir Pb-210
      d[iso] <- rbinom(1, N[iso], p[iso])
    }
    
    # Actualizar poblaciones en cadena
    N["Rn-222"] <- N["Rn-222"] - d["Rn-222"]
    N["Po-218"] <- N["Po-218"] + d["Rn-222"] - d["Po-218"]
    N["Pb-214"] <- N["Pb-214"] + d["Po-218"] - d["Pb-214"]
    N["Bi-214"] <- N["Bi-214"] + d["Pb-214"] - d["Bi-214"]
    # Po-214: instantáneo — lo que llega de Bi-214 pasa directo a Pb-210
    N["Pb-210"] <- N["Pb-210"] + d["Bi-214"]   # Po-214 no acumula
    
    trayectoria[i, ] <- N
  }
  
  suma <- suma + trayectoria
}

# Promedio punto a punto
promedio <- as.data.frame(suma / n_runs)
promedio$t <- tiempo

# Formato largo para ggplot

mc_long <- pivot_longer(promedio, -t, names_to = "isotopo", values_to = "N") |>
  # Orden de aparición en la cadena
  transform(isotopo = factor(isotopo, levels = nombres))

# Colores

colores <- c(
  "Rn-222" = "#378ADD",
  "Po-218" = "#D85A30",
  "Pb-214" = "#1D9E75",
  "Bi-214" = "#9c59d1",
  "Po-214" = "#E24B4A",
  "Pb-210" = "#BA7517"
)

# Tema y eje x compartidos
tema_base <- theme_minimal(base_size = 13) +
  theme(legend.position  = "top",
        panel.grid.minor = element_blank())

# Panel superior: Rn-222 y Pb-210 

datos_sup <- mc_long[mc_long$isotopo %in% c("Rn-222", "Pb-210"), ]

p_sup <- ggplot(datos_sup, aes(x = t, y = N, color = isotopo)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = colores) +
  scale_x_continuous(
    breaks   = seq(0, t_max, by = 25),
    sec.axis = sec_axis(~ . / 24, name = "Tiempo (días)",
                        breaks = seq(0, t_max / 24, by = 1))
  ) +
  labs(
    title    = "Cadena de desintegración del Radón-222 — Monte Carlo",
    subtitle = paste0("Rn-222 → Po-218 → Pb-214 → Bi-214 → Po-214 → Pb-210  |  ",
                      n_runs, " corridas  |  N₀ = ", N0, " átomos"),
    x     = NULL,
    y     = "Átomos N(t)",
    color = NULL
  ) +
  tema_base +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x.bottom  = element_blank(),
        axis.ticks.x.bottom = element_blank())

# Panel inferior: isótopos intermedios con escala Y propia

intermedios <- c("Po-218", "Pb-214", "Bi-214", "Po-214")
datos_inf   <- mc_long[mc_long$isotopo %in% intermedios, ]

p_inf <- ggplot(datos_inf, aes(x = t, y = N, color = isotopo)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = colores) +
  scale_x_continuous(breaks = seq(0, t_max, by = 25)) +
  labs(
    x     = "Tiempo (horas)",
    y     = "Átomos N(t)  [escala propia]",
    color = NULL
  ) +
  tema_base

# Componer con patchwork

library(patchwork)

print(p_sup / p_inf + plot_layout(heights = c(1.4, 1)))