# Codigo de validacion del metodo de Monte Carlo

# Matias Berrio
# Samuel Vargas


library(ggplot2)
library(tidyr)



T12_A  <- 60.6              # vida media Bi-212 [min]
T12_B  <- 3.05              # vida media Tl-208 [min]
lA     <- log(2) / T12_A
lB     <- log(2) / T12_B

N0     <- 1000              # átomos iniciales
t_max  <- 600               # [min]  ~10 vidas medias del padre
dt     <- 0.5               # paso de tiempo [min]
n_runs <- 30                # corridas Monte Carlo

tiempo <- seq(0, t_max, by = dt)

colores <- c("Bi-212" = "#378ADD",
             "Tl-208" = "#D85A30",
             "Pb-208" = "#1D9E75")

# 1. Solución analítica de Bateman

NA_ <- N0 * exp(-lA * tiempo)
NB_ <- N0 * lA / (lB - lA) * (exp(-lA * tiempo) - exp(-lB * tiempo))
NC_ <- N0 - NA_ - NB_

bat <- data.frame(t = tiempo, `Bi-212` = NA_, `Tl-208` = NB_, `Pb-208` = NC_,
                  check.names = FALSE) |>
  pivot_longer(-t, names_to = "isotopo", values_to = "N")

p1 <- ggplot(bat, aes(x = t, y = N, color = isotopo)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = colores) +
  labs(
    title    = "Solución analítica — Ecuaciones de Bateman",
    subtitle = "Bi-212 → Tl-208 → Pb-208 (estable)",
    x        = "Tiempo (min)", y = "Átomos N(t)", color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top",
        panel.grid.minor = element_blank())

# 2. Simulación Monte Carlo

p_A <- 1 - exp(-lA * dt)
p_B <- 1 - exp(-lB * dt)
n_steps <- length(tiempo)

set.seed(42)

# Guardar todas las corridas individualmente para Bi-212
runs_A <- matrix(0, nrow = n_steps, ncol = n_runs)
sum_A  <- numeric(n_steps)
sum_B  <- numeric(n_steps)
sum_C  <- numeric(n_steps)

for (r in seq_len(n_runs)) {
  nA <- N0; nB <- 0; nC <- 0
  vA <- numeric(n_steps); vB <- numeric(n_steps); vC <- numeric(n_steps)
  vA[1] <- nA
  
  for (i in 2:n_steps) {
    dA    <- rbinom(1, nA, p_A)
    dB    <- rbinom(1, nB, p_B)
    nA    <- nA - dA
    nB    <- nB + dA - dB
    nC    <- nC + dB
    vA[i] <- nA; vB[i] <- nB; vC[i] <- nC
  }
  runs_A[, r] <- vA
  sum_A <- sum_A + vA
  sum_B <- sum_B + vB
  sum_C <- sum_C + vC
}

# Gráfica 2: 30 trayectorias individuales del Bi-212 

corridas_df <- as.data.frame(runs_A)
colnames(corridas_df) <- paste0("r", seq_len(n_runs))
corridas_df$t <- tiempo
corridas_long <- pivot_longer(corridas_df, -t, names_to = "corrida", values_to = "N")

promedio_df <- data.frame(t = tiempo, N = rowMeans(runs_A))

p2 <- ggplot() +
  geom_line(data = corridas_long,
            aes(x = t, y = N, group = corrida),
            color = "#378ADD", alpha = 0.2, linewidth = 0.4) +
  geom_line(data = promedio_df,
            aes(x = t, y = N),
            color = "#D85A30", linewidth = 1.3) +
  annotate("text", x = 380, y = 880,
           label = "30 corridas individuales", color = "#378ADD", size = 4) +
  annotate("text", x = 380, y = 820,
           label = "Promedio punto a punto",  color = "#D85A30", size = 4) +
  labs(
    title    = "Monte Carlo — 30 trayectorias individuales del Bi-212",
    subtitle = "Cada línea azul = una historia completa de t=0 a t=600 min",
    x        = "Tiempo (min)", y = "Átomos Bi-212 N(t)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title    = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# Gráfica 3: promedio MC vs los 3 isótopos 

mc <- data.frame(
  t        = tiempo,
  `Bi-212` = sum_A / n_runs,
  `Tl-208` = sum_B / n_runs,
  `Pb-208` = sum_C / n_runs,
  check.names = FALSE
) |> pivot_longer(-t, names_to = "isotopo", values_to = "N")

p3 <- ggplot(mc, aes(x = t, y = N, color = isotopo)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = colores) +
  labs(
    title    = paste("Simulación Monte Carlo —", n_runs, "corridas (promedio)"),
    subtitle = "Bi-212 → Tl-208 → Pb-208 (estable)  |  N₀ = 1000",
    x        = "Tiempo (min)", y = "Átomos N(t)", color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title    = element_text(face = "bold"),
        legend.position = "top",
        panel.grid.minor = element_blank())

# Mostrar 

print(p1)
print(p2)
print(p3)