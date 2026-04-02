# ============================================================
# PROYECTO: ERPT EN ARGENTINA CON STVAR + LOCAL PROJECTIONS
# ============================================================

# ============================================================
# 0. LIBRERÍAS
# ============================================================

# Instalar solo si no las tenés
# install.packages(c("sstvars","ggplot2","dplyr","tseries","lpirfs"))

library(sstvars)    # STVAR
library(ggplot2)    # gráficos
library(dplyr)      # manipulación de datos
library(tseries)    # tests ADF
library(lpirfs)     # Local Projections

# ============================================================
# 1. CARGA DE DATOS
# ============================================================

# El dataset debe contener:
# fecha | CPI | TC | outputgap

datos_raw <- read.csv("datos_argentina_erpt.csv", sep = ";")

# ============================================================
# 2. TRANSFORMACIONES (CRÍTICO)
# ============================================================

# Convertimos a tasas de variación (log-diferencias)
# Esto evita problemas de no estacionariedad

datos <- datos_raw %>%
  mutate(
    infl = 100 * (log(CPI) - lag(log(CPI))),   # inflación mensual %
    depre = 100 * (log(TC) - lag(log(TC))),    # depreciación %
    outputgap = outputgap                     # asumimos ya estacionario
  ) %>%
  na.omit()  # eliminamos NA generados por lags

# ============================================================
# 3. TEST DE ESTACIONARIEDAD
# ============================================================

# Verificamos que las series sean estacionarias

print("ADF Test - Inflación")
print(adf.test(datos$infl))

print("ADF Test - Depreciación")
print(adf.test(datos$depre))

print("ADF Test - Output Gap")
print(adf.test(datos$outputgap))

# ============================================================
# 4. MATRIZ PARA STVAR
# ============================================================

# Orden IMPORTANTE:
# 1 = inflación
# 2 = depreciación
# 3 = ciclo económico

data_stvar <- as.matrix(datos[, c("infl", "depre", "outputgap")])

# ============================================================
# 5. SELECCIÓN DE LAGS (AIC)
# ============================================================

lags <- 1:4
aic_vals <- numeric(length(lags))

for (i in seq_along(lags)) {
  
  cat("Estimando modelo con lag =", lags[i], "\n")
  
  fit_tmp <- fitSTVAR(
    data = data_stvar,
    p = lags[i],
    M = 2,
    weight_function = "logistic",
    weightfun_pars = c(1,1),   # transición por inflación
    cond_dist = "Student"
  )
  
  aic_vals[i] <- AIC(fit_tmp)
}

lag_selection <- data.frame(lag = lags, AIC = aic_vals)
print(lag_selection)

# Elegir el lag con menor AIC (ej: 2)
p_opt <- lag_selection$lag[which.min(lag_selection$AIC)]

# ============================================================
# 6. ESTIMACIÓN STVAR
# ============================================================

# -----------------------------
# (A) Modelo con régimen inflacionario
# -----------------------------

fit_infl <- fitSTVAR(
  data = data_stvar,
  p = p_opt,
  M = 2,                         # 2 regímenes
  weight_function = "logistic", # transición suave
  weightfun_pars = c(1,1),      # variable 1 = inflación
  cond_dist = "Student",        # robusto a outliers
  nrounds = 10,
  seeds = 1:10,
  estim_method = "two-phase"
)

summary(fit_infl)

# -----------------------------
# (B) Modelo con régimen real (ciclo)
# -----------------------------

fit_ciclo <- fitSTVAR(
  data = data_stvar,
  p = p_opt,
  M = 2,
  weight_function = "logistic",
  weightfun_pars = c(3,1),      # variable 3 = output gap
  cond_dist = "Student",
  nrounds = 10,
  seeds = 1:10,
  estim_method = "two-phase"
)

# ============================================================
# 7. ERPT: GIRF (IMPULSE RESPONSES)
# ============================================================

# Definimos tamaño de shock realista
shock_sd <- sd(datos$depre)

# Shock en depreciación → respuesta en inflación

girf_infl <- GIRF(
  fit_infl,
  shock_size = shock_sd,
  nstep = 24,             # horizonte 24 meses
  which_shocks = 2,       # depreciación
  which_responses = 1     # inflación
)

# ============================================================
# 8. GIRF POR RÉGIMEN (CLAVE)
# ============================================================

# Régimen 1 (ej: baja inflación)
plot(girf_infl, regimes = 1)

# Régimen 2 (ej: alta inflación)
plot(girf_infl, regimes = 2)

# ============================================================
# 9. GRÁFICO PERSONALIZADO
# ============================================================

df_irf <- as.data.frame(girf_infl$irf)

ggplot(df_irf, aes(x = step, y = value)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Exchange Rate Pass-Through (ERPT)",
       subtitle = "Respuesta de inflación a shock cambiario",
       x = "Horizonte (meses)",
       y = "Respuesta de inflación (%)")

# ============================================================
# 10. DIAGNÓSTICOS
# ============================================================

# Verificación del modelo

diagnostic_plot(fit_infl)       # residuos
Portmanteau_test(fit_infl)      # autocorrelación
LR_test(fit_infl)               # lineal vs no lineal

# ============================================================
# 11. DESCOMPOSICIONES
# ============================================================

# Varianza explicada por shocks
GFEVD(fit_infl, nstep = 24)

# Descomposición histórica
hist_decomp(fit_infl)

# ============================================================
# 12. MODELO ESTRUCTURAL
# ============================================================

fit_struct <- fitSSTVAR(
  stvar = fit_infl,
  identification = "heteroskedasticity"
)

girf_struct <- GIRF(
  fit_struct,
  shock_size = shock_sd,
  nstep = 24,
  which_shocks = 2,
  which_responses = 1
)

plot(girf_struct)

# ============================================================
# 13. COMPARACIÓN: LOCAL PROJECTIONS
# ============================================================

lp_model <- lp_lin(
  endog_data = datos[, c("infl", "depre", "outputgap")],
  lags_endog_lin = p_opt,
  shock = "depre",
  response = "infl",
  hor = 24
)

plot(lp_model)

# ============================================================
# 14. INTERPRETACIÓN (GUÍA)
# ============================================================

# - Comparar magnitud de IRFs entre regímenes
# - Alta inflación → mayor ERPT esperado
# - Baja inflación → menor ERPT
# - LP = efecto promedio
# - STVAR = efecto dependiente del estado

# ============================================================
# FIN DEL SCRIPT
# ============================================================


# ============================================================
# ERPT EN ARGENTINA - STVAR + GIRF + LOCAL PROJECTIONS
# SCRIPT COMPLETO PARA TESIS / PAPER
# ============================================================

# ============================================================
# 0. LIBRERÍAS
# ============================================================

# install.packages(c("sstvars","ggplot2","dplyr","tseries","lpirfs"))

library(sstvars)
library(ggplot2)
library(dplyr)
library(tseries)
library(lpirfs)

# ============================================================
# 1. CARGA DE DATOS
# ============================================================

datos_raw <- read.csv("datos_argentina_erpt.csv", sep = ";")

# ============================================================
# 2. TRANSFORMACIONES (ESTACIONARIEDAD)
# ============================================================

datos <- datos_raw %>%
  mutate(
    infl = 100 * (log(CPI) - lag(log(CPI))),   # inflación
    depre = 100 * (log(TC) - lag(log(TC))),    # depreciación
    outputgap = outputgap
  ) %>%
  na.omit()

# ============================================================
# 3. TEST DE ESTACIONARIEDAD
# ============================================================

print(adf.test(datos$infl))
print(adf.test(datos$depre))
print(adf.test(datos$outputgap))

# ============================================================
# 4. MATRIZ STVAR
# ============================================================

data_stvar <- as.matrix(datos[, c("infl", "depre", "outputgap")])

# ============================================================
# 5. SELECCIÓN DE LAGS (AIC)
# ============================================================

lags <- 1:4
aic_vals <- numeric(length(lags))

for (i in seq_along(lags)) {
  
  fit_tmp <- fitSTVAR(
    data = data_stvar,
    p = lags[i],
    M = 2,
    weight_function = "logistic",
    weightfun_pars = c(1,1),
    cond_dist = "Student"
  )
  
  aic_vals[i] <- AIC(fit_tmp)
}

lag_selection <- data.frame(lag = lags, AIC = aic_vals)
print(lag_selection)

p_opt <- lag_selection$lag[which.min(lag_selection$AIC)]

# ============================================================
# 6. ESTIMACIÓN STVAR
# ============================================================

fit_infl <- fitSTVAR(
  data = data_stvar,
  p = p_opt,
  M = 2,
  weight_function = "logistic",
  weightfun_pars = c(1,1),
  cond_dist = "Student",
  nrounds = 10,
  seeds = 1:10,
  estim_method = "two-phase"
)

summary(fit_infl)

# ============================================================
# 7. DEFINICIÓN DEL SHOCK
# ============================================================

shock_sd <- sd(datos$depre)

# ============================================================
# 8. GIRF CON BOOTSTRAP (RESULTADO CENTRAL)
# ============================================================

nboot <- 200

girf_boot <- GIRF(
  fit_infl,
  shock_size = shock_sd,
  nstep = 24,
  which_shocks = 2,
  which_responses = 1,
  nboot = nboot,
  ci = 0.95
)

plot(girf_boot)

# ============================================================
# 9. DATAFRAME PARA GRÁFICOS
# ============================================================

df_irf <- data.frame(
  step  = 0:(length(girf_boot$irf)-1),
  mean  = girf_boot$irf,
  lower = girf_boot$lower,
  upper = girf_boot$upper
)

# ============================================================
# 10. GRÁFICO PUBLICABLE
# ============================================================

p_irf <- ggplot(df_irf, aes(x = step, y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(title = "Exchange Rate Pass-Through",
       subtitle = "Respuesta de inflación a shock cambiario",
       x = "Horizonte (meses)",
       y = "Inflación (%)")

print(p_irf)

ggsave("irf_erpt.png", plot = p_irf, width = 8, height = 5, dpi = 300)

# ============================================================
# 11. TABLA RESUMEN
# ============================================================

impacto <- df_irf$mean[1]
pico <- max(df_irf$mean)
t_pico <- which.max(df_irf$mean) - 1
persistencia <- which(abs(df_irf$mean) < 0.05)[1] - 1

tabla_resumen <- data.frame(
  Impacto = impacto,
  Pico = pico,
  Tiempo_al_pico = t_pico,
  Persistencia = persistencia
)

print(tabla_resumen)
write.csv(tabla_resumen, "tabla_erpt.csv", row.names = FALSE)

# ============================================================
# 12. GIRF POR RÉGIMEN
# ============================================================

girf_reg1 <- GIRF(
  fit_infl,
  shock_size = shock_sd,
  nstep = 24,
  which_shocks = 2,
  which_responses = 1,
  regimes = 1,
  nboot = nboot
)

girf_reg2 <- GIRF(
  fit_infl,
  shock_size = shock_sd,
  nstep = 24,
  which_shocks = 2,
  which_responses = 1,
  regimes = 2,
  nboot = nboot
)

# ============================================================
# 13. COMPARACIÓN DE REGÍMENES
# ============================================================

df_r1 <- data.frame(step = 0:(length(girf_reg1$irf)-1),
                    value = girf_reg1$irf,
                    tipo = "Baja inflación")

df_r2 <- data.frame(step = 0:(length(girf_reg2$irf)-1),
                    value = girf_reg2$irf,
                    tipo = "Alta inflación")

df_reg <- rbind(df_r1, df_r2)

ggplot(df_reg, aes(x = step, y = value, color = tipo)) +
  geom_line(size = 1.2) +
  theme_minimal() +
  labs(title = "ERPT por régimen",
       x = "Horizonte",
       y = "Inflación (%)")

# ============================================================
# 14. DIAGNÓSTICOS
# ============================================================

diagnostic_plot(fit_infl)
Portmanteau_test(fit_infl)
LR_test(fit_infl)

# ============================================================
# 15. LOCAL PROJECTIONS (COMPARACIÓN)
# ============================================================

lp_model <- lp_lin(
  endog_data = datos[, c("infl", "depre", "outputgap")],
  lags_endog_lin = p_opt,
  shock = "depre",
  response = "infl",
  hor = 24
)

plot(lp_model)

# ============================================================
# 16. ROBUSTNESS BÁSICO
# ============================================================

fit_gaussian <- fitSTVAR(
  data = data_stvar,
  p = p_opt,
  M = 2,
  weight_function = "logistic",
  weightfun_pars = c(1,1),
  cond_dist = "Gaussian"
)

# ============================================================
# FIN
# ============================================================




# ============================================================
# ERPT EN ARGENTINA - STVAR + GIRF + LOCAL PROJECTIONS
# SCRIPT CORREGIDO Y OPTIMIZADO (2026)
# ============================================================

# 0. LIBRERÍAS
# ============================================================
# install.packages(c("sstvars", "ggplot2", "dplyr", "tseries", "lpirfs"))
library(sstvars)
library(ggplot2)
library(dplyr)
library(tseries)
library(lpirfs)

# ============================================================
# 1. CARGA DE DATOS
# ============================================================
datos_raw <- read.csv("datos_argentina_erpt.csv", sep = ";")

# ============================================================
# 2. TRANSFORMACIONES (ESTACIONARIEDAD)
# ============================================================
datos <- datos_raw %>%
  mutate(
    infl  = 100 * (log(CPI) - lag(log(CPI))),      # inflación mensual %
    depre = 100 * (log(TC)  - lag(log(TC))),       # depreciación %
    outputgap = outputgap                          # ya es gap (chequeá estacionariedad)
  ) %>%
  na.omit()

# ============================================================
# 3. TEST DE ESTACIONARIEDAD (recomiendo también KPSS)
# ============================================================
adf.test(datos$infl)
adf.test(datos$depre)
adf.test(datos$outputgap)

# ============================================================
# 4. MATRIZ PARA STVAR (orden: infl, depre, outputgap)
# ============================================================
data_stvar <- as.matrix(datos[, c("infl", "depre", "outputgap")])

# ============================================================
# 5. SELECCIÓN DE LAGS (mejor usar AIC del modelo completo)
# ============================================================
lags <- 1:4
aic_vals <- numeric(length(lags))

for (i in seq_along(lags)) {
  fit_tmp <- tryCatch({
    fitSTVAR(
      data = data_stvar,
      p = lags[i],
      M = 2,
      weight_function = "logistic",
      weightfun_pars = c(1, 1),      # transición con infl lag 1
      cond_dist = "Student",
      nrounds = 8,                   # menos rondas para selección rápida
      seeds = 1:8,
      estim_method = "two-phase"
    )
  }, error = function(e) NULL)
  
  aic_vals[i] <- if (!is.null(fit_tmp)) AIC(fit_tmp) else NA
}

lag_selection <- data.frame(lag = lags, AIC = aic_vals)
print(lag_selection)
p_opt <- lag_selection$lag[which.min(lag_selection$AIC)]
cat("Lags óptimos según AIC:", p_opt, "\n")

# ============================================================
# 6. ESTIMACIÓN STVAR PRINCIPAL (Régimen inflacionario)
# ============================================================
fit_infl <- fitSTVAR(
  data = data_stvar,
  p = p_opt,
  M = 2,
  weight_function = "logistic",
  weightfun_pars = c(1, 1),          # inflación rezagada como variable de transición
  cond_dist = "Student",
  nrounds = 15,                      # aumentar si no converge bien
  seeds = 1:15,
  estim_method = "two-phase"
)

summary(fit_infl)
print(fit_infl)   # muestra parámetros clave

# ============================================================
# 7. GIRF (IMPORTANTE: corrección de sintaxis)
# ============================================================
# GIRF en sstvars funciona sobre modelos **estructurales** o reduced-form.
# Para reduced-form el shock se interpreta como generalized.

girf_obj <- GIRF(
  stvar        = fit_infl,           # objeto del modelo
  which_shocks = 2,                  # shock en variable 2 = depre
  shock_size   = 1,                  # 1 desviación estándar (recomendado)
  N            = 24,                 # horizonte (antes nstep)
  R1           = 200,                # repeticiones internas por valor inicial
  R2           = 200,                # número de valores iniciales
  init_regime  = NULL,               # NULL = mezcla de regímenes (recomendado)
  ci           = 0.95,
  ncores       = 2,                  # ajustá según tu máquina
  use_parallel = TRUE
)

plot(girf_obj)   # plot genérico del paquete

# Extraer para ggplot
df_irf <- data.frame(
  step  = 0:(length(girf_obj$girf_mean) - 1),   # chequeá el nombre exacto con str(girf_obj)
  mean  = girf_obj$girf_mean,                   # nombres pueden variar: irf, girf, etc.
  lower = girf_obj$girf_lower,
  upper = girf_obj$girf_upper
)

# ============================================================
# 8. GRÁFICO PUBLICABLE
# ============================================================
p_irf <- ggplot(df_irf, aes(x = step, y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, fill = "steelblue") +
  geom_line(color = "steelblue", linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  theme_minimal(base_size = 14) +
  labs(title = "Exchange Rate Pass-Through en Argentina",
       subtitle = "Respuesta de la inflación a un shock de depreciación (STVAR - régimen logístico)",
       x = "Horizonte (meses)",
       y = "Respuesta de Inflación (%)") +
  theme(plot.title = element_text(face = "bold"))

print(p_irf)
ggsave("ERPT_GIRF.png", p_irf, width = 9, height = 5.5, dpi = 400)

# ============================================================
# 9. TABLA RESUMEN
# ============================================================
impacto      <- df_irf$mean[2]          # mes 1 (step=0 suele ser 0)
pico         <- max(df_irf$mean, na.rm = TRUE)
t_pico       <- which.max(df_irf$mean) - 1
persistencia <- min(which(abs(df_irf$mean) < 0.05), na.rm = TRUE) - 1   # cuándo cae por debajo de 0.05pp

tabla_resumen <- data.frame(
  Impacto_inicial = impacto,
  Pico_maximo     = pico,
  Mes_del_pico    = t_pico,
  Persistencia_meses = persistencia
)
print(tabla_resumen)

# ============================================================
# 10. GIRF POR RÉGIMEN (opcional - más costoso computacionalmente)
# ============================================================
# Para ver diferencias Baja vs Alta inflación:
girf_reg1 <- GIRF(fit_infl, which_shocks = 2, shock_size = 1, N = 24,
                  init_regime = 1, R1 = 150, R2 = 150)

girf_reg2 <- GIRF(fit_infl, which_shocks = 2, shock_size = 1, N = 24,
                  init_regime = 2, R1 = 150, R2 = 150)

# ============================================================
# 11. DIAGNÓSTICOS
# ============================================================
diagnostic_plot(fit_infl)
Portmanteau_test(fit_infl)
LR_test(fit_infl)          # test linealidad vs STVAR (muy importante reportar)

# ============================================================
# 12. LOCAL PROJECTIONS (comparación lineal)
# ============================================================
# lp_lin espera data.frame, no matrix
lp_model <- lp_lin(
  endog_data     = as.data.frame(data_stvar),
  lags_endog_lin = p_opt,
  shock_type     = 1,          # 1 = unit shock (o 0 = 1 sd)
  hor            = 24,
  confint        = 1.96        # 95%
)

plot(lp_model)

# ============================================================
# 13. ROBUSTNESS (cambiar distribución y/o variable de transición)
# ============================================================
# Ejemplo: transición por ciclo económico (outputgap)
fit_ciclo <- fitSTVAR(
  data = data_stvar,
  p = p_opt,
  M = 2,
  weight_function = "logistic",
  weightfun_pars = c(3, 1),      # outputgap lag 1
  cond_dist = "Student",
  nrounds = 12,
  seeds = 1:12
)

# Gaussian como robustness
fit_gauss <- fitSTVAR(data = data_stvar, p = p_opt, M = 2,
                      weight_function = "logistic", weightfun_pars = c(1,1),
                      cond_dist = "Gaussian", nrounds = 10)

# ============================================================
# FIN


# =============================================
# GIRF GENERAL (promedio de regímenes) - ya lo tenías
# =============================================
girf_general <- GIRF(
  stvar        = fit_infl,
  which_shocks = 2,                  # shock en depreciación (variable 2)
  shock_size   = 1,                  # 1 desviación estándar (recomendado)
  N            = 24,                 # horizonte de 24 meses
  R1           = 200,                # repeticiones por valor inicial
  R2           = 200,                # cantidad de valores iniciales diferentes
  ci           = 0.95,
  ncores       = 2,                  # ajustá según tu PC
  use_parallel = TRUE
)

plot(girf_general)

# =============================================
# GIRF POR RÉGIMEN INFLACIONARIO (lo que preguntás)
# =============================================

# Régimen 1: Baja inflación (asumiendo que el régimen 1 es el de inflación baja)
girf_baja <- GIRF(
  stvar        = fit_infl,
  which_shocks = 2,
  shock_size   = 1,
  N            = 24,
  init_regime  = 1,          # ← CLAVE: condicionado a régimen 1
  R1           = 150,
  R2           = 150,
  ci           = 0.95,
  ncores       = 2,
  use_parallel = TRUE
)

# Régimen 2: Alta inflación
girf_alta <- GIRF(
  stvar        = fit_infl,
  which_shocks = 2,
  shock_size   = 1,
  N            = 24,
  init_regime  = 2,          # ← CLAVE: condicionado a régimen 2
  R1           = 150,
  R2           = 150,
  ci           = 0.95,
  ncores       = 2,
  use_parallel = TRUE
)

# Gráfico comparativo (recomendado para tesis/paper)
df_reg <- rbind(
  data.frame(step = 0:24, 
             irf = girf_baja$girf_mean,     # chequeá el nombre exacto con str(girf_baja)
             lower = girf_baja$girf_lower,
             upper = girf_baja$girf_upper,
             regimen = "Baja inflación"),
  
  data.frame(step = 0:24, 
             irf = girf_alta$girf_mean,
             lower = girf_alta$girf_lower,
             upper = girf_alta$girf_upper,
             regimen = "Alta inflación")
)

ggplot(df_reg, aes(x = step, y = irf, color = regimen, fill = regimen)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(title = "ERPT por régimen inflacionario",
       subtitle = "Respuesta de inflación a shock de depreciación (STVAR)",
       x = "Horizonte (meses)",
       y = "Respuesta (%)",
       color = "Régimen", fill = "Régimen") +
  scale_color_manual(values = c("Baja inflación" = "blue", "Alta inflación" = "red"))
# ============================================================
