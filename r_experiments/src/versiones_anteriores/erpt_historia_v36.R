cat(">>> erpt_historia v36 cargada (graficos_con_mediana: referencias por régimen con la mediana)\n")
# =====================================================================
# ERPT dependiente del estado: GIRF condicionales a la historia de cada trimestre
# Pensado para usarse después de run_vlstar() (o run_threshold_tvar()) de tu wrapper.
#
# Uso (en una celda del notebook):
#   source("erpt_condicional_historia.R")
#   out <- erpt_historia(res, data = df_ERPT_Arg_capped,
#                        endog = c("Psoja_USA", "E", "ipc", "pbird"))
#   graficar_erpt(out, archivo = "erpt_soja_p2")
# =====================================================================

library(sstvars)
library(plotly)
library(htmlwidgets)

# --- Utilidades para encontrar cosas dentro de 'res' sin conocer su estructura exacta ---
.buscar <- function(x, pred) {
  if (pred(x)) return(x)
  if (is.list(x) && !inherits(x, "stvar")) {
    for (el in x) {
      r <- .buscar(el, pred)
      if (!is.null(r)) return(r)
    }
  }
  NULL
}
.buscar_nombre <- function(x, nombre) {
  if (!is.list(x) || inherits(x, "stvar")) return(NULL)
  if (nombre %in% names(x)) return(x[[nombre]])
  for (el in x) {
    r <- .buscar_nombre(el, nombre)
    if (!is.null(r)) return(r)
  }
  NULL
}

# =====================================================================
# Cálculo
# =====================================================================
erpt_historia <- function(res,
                          data,                 # df_ERPT_Arg_capped (el mismo que usaste para estimar)
                          endog,                # mismo orden que en run_vlstar
                          shock_var  = "E",
                          resp_var   = "ipc",
                          shock_num  = NULL,    # por defecto: posición de shock_var en endog (como tu wrapper)
                          N          = 8,       # horizonte en trimestres (8 = 2 años, como tu wrapper)
                          R1         = 1000,    # simulaciones por historia
                          regimen_alto = NULL,  # régimen de mayor inflación (por defecto, el último)
                          definicion = c("impacto", "acumulado"),
                          modelo     = NULL,    # si no lo encuentra en 'res', pasalo a mano
                          escalas    = NULL,    # desvíos usados para escalar (se buscan en res$...$scaling_sd)
                          semilla    = 123,
                          ncores     = max(1, parallel::detectCores() - 1),
                          deprec_pct = NULL,     # NULL = shock de 1 desvío; un número = depreciación inicial en %
                          regimen_shock = 1) {   # régimen en el que se calibra el tamaño y el signo del shock

  definicion <- match.arg(definicion)

  # ---- Modelo estructural y escalas ----
  if (is.null(modelo)) {
    modelo <- .buscar(res, function(x) inherits(x, "stvar") &&
                        !identical(x$model$identification, "reduced_form"))
  }
  if (is.null(modelo)) stop("No encontré el modelo estructural en 'res'. Pasalo con modelo = ...")
  if (is.null(escalas)) escalas <- .buscar_nombre(res, "scaling_sd")
  if (is.null(escalas)) escalas <- sapply(data[, endog], sd, na.rm = TRUE)  # misma cuenta que fit_generic_sstvar
  M <- modelo$model$M
  if (is.null(regimen_alto)) regimen_alto <- M

  idx_e <- match(shock_var, endog)
  idx_p <- match(resp_var, endog)
  if (is.null(shock_num)) shock_num <- idx_e
  sd_e <- unname(escalas[idx_e])
  sd_p <- unname(escalas[idx_p])

  Y     <- as.matrix(modelo$data)
  p     <- modelo$model$p
  d     <- ncol(Y)
  T_obs <- nrow(Y)
  if (T_obs != nrow(data)) {
    stop("El modelo tiene ", T_obs, " observaciones y 'data' tiene ", nrow(data),
         ". Revisá que sea el mismo data.frame (sin filas eliminadas).")
  }

  # ---- Fechas trimestrales ----
  q        <- as.integer(substr(data$trimestre, 1, 1))
  anio_dec <- data$año + (q - 1) / 4
  etiqueta <- paste0(q, "q", data$año)

  # ---- Una historia de longitud p por trimestre; el shock ocurre en t+1 ----
  t_fin <- p:(T_obs - 1)
  K     <- length(t_fin)
  init_arr <- array(NA_real_, dim = c(p, d, K))
  for (k in seq_len(K)) {
    init_arr[, , k] <- Y[(t_fin[k] - p + 1):t_fin[k], , drop = FALSE]  # última fila = dato más reciente
  }
  t_shock <- t_fin + 1

  # Peso del régimen alto en el trimestre del shock
  # (fila j de transition_weights corresponde a la observación p + j)
  pesos   <- modelo$transition_weights[t_shock - p, , drop = FALSE]
  alpha   <- pesos[, regimen_alto]
  regimen <- apply(pesos, 1, which.max)

  # ---- GIRF condicional a cada historia ----
  # Mismo shock (1 desvío estructural) y misma semilla en todas las historias:
  # la variación entre trimestres refleja solo el estado de la economía.
  tam <- .tam_shock(modelo, escalas, idx_e, shock_num, deprec_pct, regimen_shock)
  message("Estimando ", K, " GIRFs (una por trimestre), ", R1, " simulaciones cada una...")
  g <- GIRF(modelo,
            which_shocks     = shock_num,
            shock_size       = tam,
            N                = N,
            R1               = R1,
            R2               = K,            # debe coincidir con dim(init_values)[3]
            init_values      = init_arr,
            which_cumulative = numeric(0),   # acumulamos a mano abajo
            seeds            = rep(semilla, K),
            ncores           = ncores)

  resp <- g$all_girfs[[1]]               # [horizonte, variable, historia]
  if (is.null(resp) || dim(resp)[3] != K) stop("GIRF no devolvió una respuesta por historia.")

  r_e <- matrix(resp[, idx_e, ], nrow = N + 1)   # unidades escaladas
  r_p <- matrix(resp[, idx_p, ], nrow = N + 1)
  cum_p <- apply(r_p, 2, cumsum)
  cum_e <- apply(r_e, 2, cumsum)

  denom_esc <- if (definicion == "impacto") matrix(r_e[1, ], N + 1, K, byrow = TRUE) else cum_e

  erpt_escalado <- cum_p / denom_esc                  # comparable con tu print_passthrough actual
  erpt          <- erpt_escalado * (sd_p / sd_e)      # en unidades originales (lo que va a la tesis)

  dimnames(erpt) <- dimnames(erpt_escalado) <-
    list(horizonte = 0:N, trimestre = etiqueta[t_shock])

  # ---- Chequeo: promedio del ERPT en las historias dominadas por cada régimen ----
  chequeo <- .chequeo_regimenes(erpt, pesos, N)
  cat("\nTrimestres dominados por cada régimen (peso > 0.9):",
      paste0("régimen ", seq_len(M), " = ", colSums(pesos > 0.9), collapse = " | "), "\n")
  cat("Shock usado:", shock_num, "| Factor de conversión a unidades originales (sd_", resp_var,
      " / sd_", shock_var, "): ", round(sd_p / sd_e, 4), "\n\n", sep = "")
  cat("Promedio entre trimestres:\n")
  print(round(chequeo, 4))
  chequeo_mediana <- .chequeo_mediana(erpt, pesos, N)
  cat("\nMediana entre trimestres (robusta a trimestres con respuestas extremas):\n")
  print(round(chequeo_mediana, 4))

  invisible(list(erpt = erpt, erpt_escalado = erpt_escalado, alpha = alpha,
                 chequeo_mediana = chequeo_mediana,
                 pesos = pesos, regimen = regimen, M = M, shock_num = shock_num,
                 anio_dec = anio_dec[t_shock], etiqueta = etiqueta[t_shock],
                 N = N, definicion = definicion, chequeo = chequeo,
                 factor = sd_p / sd_e))
}


# --- Ayudas para trabajar con 2 o 3 regímenes ---
# Completa objetos out_... creados con versiones anteriores (solo 2 regímenes)
.normalizar_out <- function(out) {
  if (is.null(out$pesos)) {
    out$pesos   <- cbind(1 - out$alpha, out$alpha)
    out$M       <- 2
    out$regimen <- ifelse(out$alpha > 0.5, 2, 1)
  }
  out
}
# Promedio del ERPT en las historias dominadas por cada régimen (peso > 0.9)
# Mediana del ERPT entre las historias dominadas por cada régimen (robusta a historias extremas)
.chequeo_mediana <- function(mat, pesos, N) {
  M <- ncol(pesos)
  df <- data.frame(horizonte = 0:N)
  nombres <- if (M == 2) c("reg_bajo_mediana", "reg_alto_mediana") else paste0("reg_", seq_len(M), "_mediana")
  for (m in seq_len(M)) {
    sel <- pesos[, m] > 0.9
    df[[nombres[m]]] <- if (any(sel)) apply(mat[, sel, drop = FALSE], 1, median) else NA
  }
  df
}
.chequeo_regimenes <- function(mat, pesos, N) {
  M <- ncol(pesos)
  df <- data.frame(horizonte = 0:N)
  nombres <- if (M == 2) c("reg_bajo_original", "reg_alto_original") else paste0("reg_", seq_len(M), "_original")
  for (m in seq_len(M)) df[[nombres[m]]] <- rowMeans(mat[, pesos[, m] > 0.9, drop = FALSE])
  df
}
# Estadístico de referencia de los gráficos: "promedio" (por defecto) o "mediana"
.estadistico <- function(out) if (identical(out$estadistico, "mediana")) "mediana" else "promedio"
.ref_por_regimen <- function(out, mat = out$erpt, N = out$N) {
  if (.estadistico(out) == "mediana") .chequeo_mediana(mat, out$pesos, N) else .chequeo_regimenes(mat, out$pesos, N)
}
.nombre_est <- function(out, mayus = FALSE) {
  e <- .estadistico(out); if (mayus) paste0(toupper(substr(e, 1, 1)), substr(e, 2, nchar(e))) else e
}

# ---------------------------------------------------------------------
# graficos_con_mediana(out): las líneas y tablas de referencia por régimen de
# todos los gráficos pasan a ser la MEDIANA de los trimestres de cada régimen
# (peso > 0,9) en lugar del promedio. Usar en los modelos que se reportan con
# la mediana (TVAR-indt, TVAR 3 reg. p=3 sin cappeo). Se puede combinar con
# acotar_para_graficar() en cualquier orden: acotar no cambia la mediana.
# ---------------------------------------------------------------------
graficos_con_mediana <- function(out) {
  out$estadistico <- "mediana"
  cat("Gráficos con la MEDIANA por régimen como referencia.\n")
  out
}
# Colores por régimen
.colores_reg <- function(M) {
  if (M == 2) c("steelblue", "firebrick")
  else if (M == 3) c("steelblue", "orange", "firebrick")
  else grDevices::hcl.colors(M, "Temps")
}
# Valor entre 0 y 1 para colorear cada historia: alfa (2 regímenes) o régimen (3 o más)
.valor_color <- function(out) {
  out <- .normalizar_out(out)
  if (out$M == 2) out$alpha else (out$regimen - 1) / (out$M - 1)
}
.escala_plotly <- function(M) {
  if (M == 2) return(list(c(0, "steelblue"), c(1, "firebrick")))
  cols <- .colores_reg(M); esc <- list()
  for (m in seq_len(M)) {
    esc[[length(esc) + 1]] <- c((m - 1) / M, cols[m])
    esc[[length(esc) + 1]] <- c(m / M - 1e-6, cols[m])
  }
  esc[[length(esc)]][1] <- 1
  esc
}
.barra_plotly <- function(M) {
  if (M == 2) return(list(title = "\u03b1 régimen alto"))
  list(title = "Régimen", tickvals = (seq_len(M) - 1) / (M - 1),
       ticktext = paste("Régimen", seq_len(M)))
}
.color_historia <- function(out, k) {
  v <- .valor_color(out)[k]
  cols <- .colores_reg(.normalizar_out(out)$M)
  rgb(colorRamp(cols)(v), maxColorValue = 255)
}
.texto_estado <- function(out, k) {
  out <- .normalizar_out(out)
  if (out$M == 2) paste0("alfa régimen alto = ", sprintf("%.2f", out$alpha[k]))
  else paste0("régimen ", out$regimen[k])
}

# Tipo de función de transición ("threshold" o "logistic")
.tipo_transicion <- function(modelo) {
  tipo <- modelo$model$weight_function
  if (is.null(tipo)) {
    tipo <- if (all(modelo$transition_weights %in% c(0, 1))) "threshold" else "logistic"
  }
  tipo
}

# --- Distribución de los errores del modelo ("Gaussian", "Student", ...) ---
.dist_errores <- function(modelo) {
  cd <- modelo$model$cond_dist
  if (is.null(cd)) "Gaussian" else cd
}

# --- Parámetros de la distribución (grados de libertad) al final de params ---
.n_param_dist <- function(modelo) {
  d <- ncol(as.matrix(modelo$data))
  switch(.dist_errores(modelo), Student = 1, ind_Student = d, ind_skewed_t = 2 * d, 0)
}

# --- Parámetros de la función de transición, sin los de la distribución ---
# Orden de params en sstvars: (..., parámetros de transición, parámetros de la distribución)
.param_transicion <- function(modelo) {
  pr <- modelo$params
  nd <- .n_param_dist(modelo)
  if (nd > 0) pr <- head(pr, length(pr) - nd)
  M <- modelo$model$M
  n_w <- switch(.tipo_transicion(modelo),
                threshold   = M - 1,
                logistic    = 2,
                exponential = 2,
                0)
  if (identical(.tipo_transicion(modelo), "mlogit")) {
    wp <- modelo$model$weightfun_pars
    k  <- if (is.list(wp) && !is.null(wp$vars)) 1 + length(wp$vars) * wp$lags else 2
    n_w <- (M - 1) * k
  }
  tail(pr, n_w)
}

# --- Grados de libertad de la t de Student (NA si el modelo es gaussiano) ---
.grados_libertad <- function(modelo) {
  if (.dist_errores(modelo) != "Student") return(NA)
  tail(modelo$params, 1)
}

# --- Cambios de régimen dominante según el valor de la variable de transición ---
# Para modelos sin umbral explícito (mlogit): ordena los trimestres por la variable
# de transición (z, unidades originales) y marca dónde cambia el régimen de mayor peso.
.cruces_empiricos <- function(z, pesos) {
  o   <- order(z)
  dom <- apply(pesos[o, , drop = FALSE], 1, which.max)
  zs  <- z[o]
  k   <- which(diff(dom) != 0)
  if (!length(k)) return(data.frame(de = integer(0), a = integer(0), z = numeric(0)))
  data.frame(de = dom[k], a = dom[k + 1], z = (zs[k] + zs[k + 1]) / 2)
}


# --- Tamaño del shock para las GIRF ---
# El signo de un shock identificado es arbitrario: se elige el signo que DEPRECIA E en el
# impacto (régimen 1), para medir siempre el ERPT de una depreciación (en un modelo no lineal,
# una apreciación no tiene por qué dar lo mismo). Con deprec_pct se elige además el tamaño:
# el shock que deprecia E en deprec_pct % en el impacto del régimen 1.
.tam_shock <- function(modelo, escalas, idx_e, shock_num, deprec_pct = NULL, regimen_shock = 1) {
  imp_e <- linear_IRF(modelo, N = 1, regime = regimen_shock)$point_est[idx_e, shock_num, 1] * unname(escalas[idx_e])
  tam <- if (is.null(deprec_pct)) 1 else log(1 + deprec_pct / 100) / abs(imp_e)
  tam <- tam * sign(imp_e)
  message(sprintf("Shock: %.3f desvíos estructurales%s; en el impacto (régimen %d) deprecia E en %.1f%%.",
                  abs(tam), if (tam < 0) " con el signo invertido" else "",
                  regimen_shock, 100 * (exp(abs(tam * imp_e)) - 1)))
  if (.dist_errores(modelo) %in% c("ind_Student", "ind_skewed_t") && modelo$model$M > 1) {
    for (m in seq_len(modelo$model$M)) {
      im <- linear_IRF(modelo, N = 1, regime = m)$point_est[idx_e, shock_num, 1] * unname(escalas[idx_e])
      message(sprintf("   régimen %d: el mismo shock mueve E en %+.1f%% en el impacto", m, 100 * (exp(tam * im) - 1)))
    }
  }
  tam
}

# --- PDF con soporte para la letra griega alfa (cairo si está disponible) ---
.abrir_pdf <- function(archivo, width, height) {
  if (isTRUE(capabilities("cairo"))) cairo_pdf(archivo, width = width, height = height)
  else pdf(archivo, width = width, height = height)
}
.letra_alfa <- function() if (isTRUE(capabilities("cairo"))) "\u03b1" else "alfa"

# --- Etiquetas trimestrales tipo "1q2002" a partir del año decimal ---
.etiqueta_trim <- function(x) {
  a <- floor(x + 1e-9); q <- round((x - a) * 4) + 1
  paste0(q, "q", a)
}
.ticks_trim <- function(x, cada = 10) {
  a <- floor(x + 1e-9); q <- round((x - a) * 4) + 1
  sel <- q == 1 & a %% cada == 0
  list(at = x[sel], lab = .etiqueta_trim(x[sel]))
}

# =====================================================================
# Gráficos
# =====================================================================
graficar_erpt <- function(out, archivo = "erpt", horizontes_2d = c(4, out$N),
                          titulo = "ERPT dependiente del estado", eje_z = "ERPT") {
  out <- .normalizar_out(out)
  N <- out$N
  M <- out$M
  write.csv(out$erpt, paste0(archivo, "_matriz.csv"))
  write.csv(data.frame(trimestre = out$etiqueta, regimen = out$regimen,
                       round(out$pesos, 4)),
            paste0(archivo, "_regimenes.csv"), row.names = FALSE)

  # ---- 3D interactivo (color = alfa con 2 regímenes; régimen con 3 o más) ----
  color_alpha <- matrix(rep(.valor_color(out), each = N + 1), nrow = N + 1)
  tk <- .ticks_trim(out$anio_dec)
  # Texto completo de cada punto (plotly no reemplaza %{text} en superficies 3D)
  etiq_mat    <- matrix(rep(.etiqueta_trim(out$anio_dec), each = N + 1), nrow = N + 1)
  h_mat       <- matrix(rep(0:N, times = ncol(out$erpt)), nrow = N + 1)
  texto_hover <- matrix(paste0(etiq_mat, "<br>Horizonte: ", h_mat,
                               "<br>", eje_z, ": ", sprintf("%.3f", out$erpt)),
                        nrow = N + 1)
  g3d <- plot_ly(x = out$anio_dec, y = 0:N, z = out$erpt, type = "surface",
                 text = texto_hover,
                 hoverinfo = "text",
                 surfacecolor = color_alpha, cmin = 0, cmax = 1,
                 colorscale = .escala_plotly(M),
                 colorbar = .barra_plotly(M)) %>%
    layout(title = titulo,
           scene = list(xaxis = list(title = "Trimestre del shock",
                                     tickvals = tk$at, ticktext = tk$lab),
                        yaxis = list(title = "Horizonte (trimestres)"),
                        zaxis = list(title = eje_z)))

  # ---- 3D estático para la tesis ----
  .abrir_pdf(paste0(archivo, "_3d.pdf"), width = 9, height = 7)
  persp(x = out$anio_dec, y = 0:N, z = t(out$erpt),
        theta = -40, phi = 25, expand = 0.6, ticktype = "detailed",
        col = "lightsteelblue", border = NA, shade = 0.4,
        xlab = "Trimestre del shock", ylab = "Horizonte (trimestres)", zlab = eje_z,
        main = titulo)
  dev.off()

  # ---- 2D complementario: ERPT a horizontes fijos + alpha sombreado ----
  dibujar_2d <- function() {
    par(mar = c(4, 4, 3, 4))
    rango <- range(out$erpt[horizontes_2d + 1, ], na.rm = TRUE)
    plot(NA, xlim = range(out$anio_dec), ylim = rango, xaxt = "n",
         xlab = "Trimestre del shock", ylab = eje_z, main = paste(titulo, "- horizontes seleccionados"))
    axis(1, at = tk$at, labels = tk$lab)
    vc <- .valor_color(out)
    polygon(c(out$anio_dec, rev(out$anio_dec)),
            c(rango[1] + vc * diff(rango), rep(rango[1], length(vc))),
            col = adjustcolor("firebrick", 0.15), border = NA)
    cols <- c("steelblue", "black", "darkgreen", "purple")
    for (i in seq_along(horizontes_2d)) {
      lines(out$anio_dec, out$erpt[horizontes_2d[i] + 1, ], col = cols[i], lwd = 1.5)
    }
    if (M == 2) {
      axis(4, at = rango[1] + c(0, 0.5, 1) * diff(rango), labels = c(0, 0.5, 1))
      mtext(paste(.letra_alfa(), "régimen alto (sombreado)"), side = 4, line = 2.5)
    } else {
      axis(4, at = rango[1] + (seq_len(M) - 1) / (M - 1) * diff(rango), labels = seq_len(M))
      mtext("Régimen (sombreado)", side = 4, line = 2.5)
    }
    legend("topleft", bty = "n", lwd = 1.5, col = cols[seq_along(horizontes_2d)],
           legend = paste0("h = ", horizontes_2d, " trimestres"))
  }
  dibujar_2d()                                   # en el notebook
  .abrir_pdf(paste0(archivo, "_2d.pdf"), width = 10, height = 5)
  dibujar_2d()                                   # en el PDF
  dev.off()

  # ---- Guardar el 3D interactivo (al final, para que un problema acá no impida los PDF) ----
  tryCatch(
    saveWidget(g3d, paste0(archivo, "_3d.html"), selfcontained = FALSE),
    error = function(e) message("No se pudo guardar el HTML: ", conditionMessage(e))
  )

  g3d
}


# =====================================================================
# Variable de transición, umbral y pesos alpha a lo largo del tiempo
# =====================================================================
# Panel superior: inflación del trimestre ANTERIOR (la variable de transición,
# en unidades originales), con el umbral c y la franja donde alpha va de 0,1 a 0,9.
# Panel inferior: alpha del régimen alto en cada trimestre.
# Como alpha_t depende de la inflación de t-1, ambos paneles quedan alineados:
# cuando la línea cruza el umbral, alpha pasa por 0,5.
graficar_transicion <- function(res, data, endog, trans_var = "ipc",
                                regimen_alto = NULL, archivo = "transicion",
                                c_escalado = NULL, gamma = NULL,
                                modelo = NULL, escalas = NULL) {
  if (is.null(modelo)) {
    modelo <- .buscar(res, function(x) inherits(x, "stvar") &&
                        !identical(x$model$identification, "reduced_form"))
  }
  if (is.null(modelo)) stop("No encontré el modelo estructural en 'res'. Pasalo con modelo = ...")
  if (is.null(escalas)) escalas <- .buscar_nombre(res, "scaling_sd")
  if (is.null(escalas)) escalas <- sapply(data[, endog], sd, na.rm = TRUE)

  M     <- modelo$model$M
  if (is.null(regimen_alto)) regimen_alto <- M
  tipo  <- .tipo_transicion(modelo)
  es_umbral <- identical(tipo, "threshold")
  es_mlogit <- identical(tipo, "mlogit")
  idx   <- match(trans_var, endog)
  sd_z  <- unname(escalas[idx])
  Y     <- as.matrix(modelo$data)
  p     <- modelo$model$p
  T_obs <- nrow(Y)

  q        <- as.integer(substr(data$trimestre, 1, 1))
  anio_dec <- data$año + (q - 1) / 4
  t_idx <- (p + 1):T_obs
  z     <- Y[t_idx - 1, idx] * sd_z                  # inflación de t-1, unidades originales
  pesos <- modelo$transition_weights[t_idx - p, , drop = FALSE]
  fecha <- anio_dec[t_idx]

  if (es_umbral) {
    # En un TVAR gaussiano, los últimos M-1 parámetros son los umbrales
    umbrales <- sort(.param_transicion(modelo)) * sd_z
    regimen  <- apply(pesos, 1, which.max)
    cat("Umbrales (inflación trimestral de t-1):",
        paste0(round(100 * (exp(umbrales) - 1), 1), "% trimestral (",
               round(100 * (exp(4 * umbrales) - 1), 0), "% anualizada)", collapse = " | "), "\n")
    cat("Trimestres por régimen:", paste0("régimen ", seq_len(M), " = ", tabulate(regimen, M),
                                          collapse = " | "), "\n")
  } else if (es_mlogit) {
    regimen <- apply(pesos, 1, which.max)
    cr <- .cruces_empiricos(z, pesos)
    umbrales <- cr$z
    cat("Transición mlogit: cambios de régimen dominante (inflación trimestral de t-1):\n")
    if (nrow(cr)) print(data.frame(de_regimen = cr$de, a_regimen = cr$a,
                                   inflacion_trim_pct = round(100 * (exp(cr$z) - 1), 1)),
                        row.names = FALSE)
    cat("Trimestres por régimen (mayor peso):", paste0("régimen ", seq_len(M), " = ",
        tabulate(regimen, M), collapse = " | "), "\n")
  } else {
    if (is.null(c_escalado) || is.null(gamma)) {
      wp <- .param_transicion(modelo)
      if (is.null(c_escalado)) c_escalado <- wp[1]
      if (is.null(gamma))      gamma      <- wp[2]
    }
    cat("Parámetros usados -> c (escalado):", round(c_escalado, 2),
        "| gamma:", round(gamma, 2),
        "\n(deben coincidir con 'Weight params' de print(modelo); si no, pasalos a mano)\n")
    c_orig <- c_escalado * sd_z
    ancho  <- log(9) / gamma * sd_z
    alpha  <- pesos[, regimen_alto]
    cat("Umbral en unidades originales:", round(c_orig, 4),
        "| franja 10%-90%:", round(c_orig - ancho, 4), "a", round(c_orig + ancho, 4), "\n")
  }

  dibujar <- function() {
    op <- par(mfrow = c(2, 1), mar = c(2, 4.5, 2.5, 1), oma = c(2, 0, 0, 0))
    on.exit(par(op))
    tk <- .ticks_trim(fecha)
    plot(fecha, z, type = "n", xlab = "", ylab = "Inflación trimestral (t-1)", xaxt = "n",
         main = if (es_umbral) "Variable de transición y umbrales estimados"
                else if (es_mlogit) "Variable de transición y cambios de régimen dominante"
                else "Variable de transición y umbral estimado")
    axis(1, at = tk$at, labels = tk$lab)
    if (es_umbral || es_mlogit) {
      lines(fecha, z, col = "black", lwd = 1.2)
      abline(h = umbrales, col = "firebrick", lty = 2, lwd = 1.5)
      legend("topleft", bty = "n", cex = 0.8,
             legend = c("Inflación (variable de transición)",
                        if (es_mlogit) "Cambio de régimen dominante" else "Umbrales"),
             col = c("black", "firebrick"), lty = c(1, 2), lwd = c(1.2, 1.5))
    } else {
      rect(min(fecha), c_orig - ancho, max(fecha), c_orig + ancho,
           col = adjustcolor("grey60", 0.3), border = NA)
      lines(fecha, z, col = "black", lwd = 1.2)
      abline(h = c_orig, col = "firebrick", lty = 2, lwd = 1.5)
      legend("topleft", bty = "n", cex = 0.8,
             legend = c("Inflación (variable de transición)", "Umbral c",
                        paste("Franja", .letra_alfa(), "entre 0,1 y 0,9")),
             col = c("black", "firebrick", adjustcolor("grey60", 0.6)),
             lty = c(1, 2, NA), lwd = c(1.2, 1.5, NA), pch = c(NA, NA, 15), pt.cex = 2)
    }

    if (es_umbral) {
      cols <- .colores_reg(M)
      plot(fecha, regimen, type = "n", ylim = c(0.5, M + 0.5), yaxt = "n", xaxt = "n", xlab = "",
           ylab = "Régimen", main = "Régimen en cada trimestre")
      axis(1, at = tk$at, labels = tk$lab)
      axis(2, at = seq_len(M), labels = seq_len(M), las = 1)
      segments(fecha, 0.5, fecha, regimen, col = adjustcolor(cols[regimen], 0.6), lwd = 2)
      lines(fecha, regimen, type = "s", col = "grey30")
    } else if (es_mlogit) {
      cols <- .colores_reg(M)
      plot(fecha, pesos[, 1], type = "n", ylim = c(0, 1), xlab = "", xaxt = "n",
           ylab = "Peso", main = "Peso de cada régimen")
      axis(1, at = tk$at, labels = tk$lab)
      for (m in seq_len(M)) lines(fecha, pesos[, m], col = cols[m], lwd = 1.5)
      legend("topleft", bty = "n", cex = 0.8, legend = paste("Régimen", seq_len(M)),
             col = cols, lwd = 1.5)
    } else {
      plot(fecha, alpha, type = "n", ylim = c(0, 1), xlab = "", xaxt = "n",
           ylab = paste(.letra_alfa(), "régimen alto"), main = "Peso del régimen de alta inflación")
      axis(1, at = tk$at, labels = tk$lab)
      polygon(c(fecha, rev(fecha)), c(alpha, rep(0, length(alpha))),
              col = adjustcolor("firebrick", 0.35), border = "firebrick")
      abline(h = c(0.1, 0.9), col = "grey40", lty = 3)
    }
    mtext("Trimestre", side = 1, outer = TRUE, line = 0.5)
  }

  dibujar()                                        # en el notebook
  .abrir_pdf(paste0(archivo, ".pdf"), width = 10, height = 7)
  dibujar()                                        # en el PDF
  dev.off()

  invisible(data.frame(trimestre = .etiqueta_trim(fecha), fecha = fecha,
                       inflacion_t_menos_1 = z, round(pesos, 4)))
}


# =====================================================================
# ERPT incremental: cuánto aumenta el ERPT en cada trimestre
# (ERPT_h - ERPT_{h-1}; en h = 0 es el propio ERPT de impacto)
# No recalcula GIRFs: usa el resultado de erpt_historia() que ya tenés.
# =====================================================================
erpt_incremental <- function(out) {
  inc <- rbind(out$erpt[1, , drop = FALSE],
               out$erpt[-1, , drop = FALSE] - out$erpt[-nrow(out$erpt), , drop = FALSE])
  dimnames(inc) <- dimnames(out$erpt)
  inc_esc <- rbind(out$erpt_escalado[1, , drop = FALSE],
                   out$erpt_escalado[-1, , drop = FALSE] -
                     out$erpt_escalado[-nrow(out$erpt_escalado), , drop = FALSE])
  dimnames(inc_esc) <- dimnames(out$erpt_escalado)

  out <- .normalizar_out(out)
  res_inc <- out
  res_inc$erpt <- inc
  res_inc$erpt_escalado <- inc_esc
  res_inc$chequeo <- .chequeo_regimenes(inc, out$pesos, out$N)
  print(round(res_inc$chequeo, 4))
  invisible(res_inc)
}


# =====================================================================
# GIRF por régimen en UNIDADES ORIGINALES
# Reemplaza a los gráficos de GIRF de sstvars (que están en desvíos estándar).
# Cada respuesta se expresa por cada 1% de depreciación inicial de E:
#   - panel de E: depreciación acumulada (vale 1 en h = 0)
#   - panel de ipc: ERPT acumulado relativo al shock inicial
# Bandas: por defecto 68% y 95% entre historias iniciales del régimen
# (misma lógica que las bandas de sstvars; no incluyen incertidumbre de parámetros).
# =====================================================================
graficar_girf_original <- function(res, endog, shock_var = "E", shock_num = NULL,
                                   N = 8, R1 = 500, R2 = 250, regimenes = 1:2,
                                   acumular = TRUE, archivo = "girf_original",
                                   bandas = c(0.68, 0.95), variables = c("ipc", "pbird"),
                                   modelo = NULL, escalas = NULL, semilla = 1,
                                   ncores = max(1, parallel::detectCores() - 1),
                                   deprec_pct = NULL, regimen_shock = 1) {
  if (is.null(modelo)) {
    modelo <- .buscar(res, function(x) inherits(x, "stvar") &&
                        !identical(x$model$identification, "reduced_form"))
  }
  if (is.null(modelo)) stop("No encontré el modelo estructural en 'res'. Pasalo con modelo = ...")
  if (is.null(escalas)) escalas <- .buscar_nombre(res, "scaling_sd")
  if (is.null(escalas)) stop("No encontré los desvíos de escalado en 'res'. Pasalos con escalas = ...")

  d     <- length(endog)
  idx_e <- match(shock_var, endog)
  if (is.null(shock_num)) shock_num <- idx_e

  tam <- .tam_shock(modelo, escalas, idx_e, shock_num, deprec_pct, regimen_shock)
  resultados <- list()
  for (r in regimenes) {
    message("GIRF del régimen ", r, "...")
    g <- GIRF(modelo, which_shocks = shock_num, shock_size = tam, N = N,
              R1 = R1, R2 = R2, init_regime = r, which_cumulative = numeric(0),
              seeds = semilla + seq_len(R2), ncores = ncores)
    a <- g$all_girfs[[1]][, seq_len(d), , drop = FALSE]      # [N+1, d, R2], escalado
    for (j in seq_len(d)) a[, j, ] <- a[, j, ] * unname(escalas[j])   # -> unidades originales
    if (acumular) a <- apply(a, c(2, 3), cumsum)              # respuestas acumuladas
    a <- sweep(a, 3, a[1, idx_e, ], "/")                      # por 1% de depreciación inicial
    resultados[[r]] <- a
  }

  cols <- c("steelblue", "firebrick", "darkgreen", "purple")
  resumen <- list()
  sel <- match(variables, endog)
  if (anyNA(sel)) stop("Alguna de 'variables' no está en 'endog'.")
  nfil <- ceiling(length(sel) / 2); ncol_g <- min(2, length(sel))
  dibujar <- function() {
    op <- par(mfrow = c(nfil, ncol_g), mar = c(4, 4.5, 2.5, 1))
    on.exit(par(op))
    for (j in sel) {
      bnd <- lapply(regimenes, function(r) {
        m <- resultados[[r]][, j, ]
        bi <- min(bandas); be <- max(bandas)
        list(med = rowMeans(m),
             q05 = apply(m, 1, quantile, (1 - be) / 2), q95 = apply(m, 1, quantile, 1 - (1 - be) / 2),
             q16 = apply(m, 1, quantile, (1 - bi) / 2), q84 = apply(m, 1, quantile, 1 - (1 - bi) / 2))
      })
      rango <- range(unlist(lapply(bnd, function(b) c(b$q05, b$q95))))
      plot(0:N, rep(0, N + 1), type = "n", ylim = rango,
           xlab = "Horizonte (trimestres)",
           ylab = if (acumular) "Resp. acumulada por 1% de dep." else "Respuesta por 1% de dep.",
           main = endog[j])
      abline(h = 0, col = "grey60")
      for (k in seq_along(regimenes)) {
        b <- bnd[[k]]
        polygon(c(0:N, N:0), c(b$q05, rev(b$q95)), col = adjustcolor(cols[k], 0.12), border = NA)
        polygon(c(0:N, N:0), c(b$q16, rev(b$q84)), col = adjustcolor(cols[k], 0.25), border = NA)
        lines(0:N, b$med, col = cols[k], lwd = 2)
      }
      if (j == sel[1]) legend("topleft", bty = "n",
                         lwd = c(rep(2, length(regimenes)), NA),
                         col = c(cols[seq_along(regimenes)], NA),
                         legend = c(paste("Régimen", regimenes),
                                    paste0("Bandas: ", round(100 * min(bandas)), "% y ",
                                           round(100 * max(bandas)), "%")))
    }
  }

  dibujar()                                            # en el notebook
  .abrir_pdf(paste0(archivo, ".pdf"), width = 5 * ncol_g, height = 4 * nfil)
  dibujar()                                            # en el PDF
  dev.off()

  # Tabla: promedio por régimen, unidades originales, por 1% de depreciación inicial
  for (r in regimenes) {
    tab <- sapply(sel, function(j) rowMeans(resultados[[r]][, j, , drop = FALSE][, 1, ]))
    tab <- matrix(tab, nrow = N + 1)
    colnames(tab) <- endog[sel]
    tab <- data.frame(horizonte = 0:N, round(tab, 4), check.names = FALSE)
    cat("\n+++ Régimen", r, "- respuestas", if (acumular) "acumuladas" else "",
        "por 1% de depreciación inicial (unidades originales) +++\n")
    print(tab, row.names = FALSE)
    resumen[[r]] <- tab
  }
  invisible(list(respuestas = resultados, resumen = resumen))
}


# =====================================================================
# IRF LINEAL por régimen en UNIDADES ORIGINALES (los dos regímenes juntos)
# Respuestas por cada 1% de depreciación inicial de E.
# Nota: sstvars no ofrece bandas de confianza para IRF lineales cuando la
# dinámica autorregresiva cambia entre regímenes (tu caso), así que se
# grafican solo las estimaciones puntuales.
# =====================================================================
graficar_irf_original <- function(res, endog, shock_var = "E", shock_num = NULL,
                                  N = 8, regimenes = 1:2, acumular = TRUE,
                                  archivo = "irf_original", modelo = NULL, escalas = NULL,
                                  variables = c("ipc", "pbird")) {
  if (is.null(modelo)) {
    modelo <- .buscar(res, function(x) inherits(x, "stvar") &&
                        !identical(x$model$identification, "reduced_form"))
  }
  if (is.null(modelo)) stop("No encontré el modelo estructural en 'res'. Pasalo con modelo = ...")
  if (is.null(escalas)) escalas <- .buscar_nombre(res, "scaling_sd")
  if (is.null(escalas)) stop("No encontré los desvíos de escalado en 'res'. Pasalos con escalas = ...")

  d     <- length(endog)
  idx_e <- match(shock_var, endog)
  if (is.null(shock_num)) shock_num <- idx_e

  resultados <- list()
  for (r in regimenes) {
    irf <- linear_IRF(modelo, N = N, regime = r, which_cumulative = numeric(0))
    pe  <- irf$point_est
    dims <- dim(pe)
    hdim <- which(dims == N + 1)
    if (length(hdim) != 1) stop("No pude identificar la dimensión del horizonte en point_est.")
    pe <- aperm(pe, c(hdim, setdiff(1:3, hdim)))       # [horizonte, variable, shock]
    m  <- pe[, seq_len(d), shock_num]                   # [N+1, d], escalado
    for (j in seq_len(d)) m[, j] <- m[, j] * unname(escalas[j])   # -> unidades originales
    if (acumular) m <- apply(m, 2, cumsum)
    m <- m / m[1, idx_e]                                # por 1% de depreciación inicial
    colnames(m) <- endog
    resultados[[r]] <- m
  }

  cols <- c("steelblue", "firebrick", "darkgreen", "purple")
  sel <- match(variables, endog)
  if (anyNA(sel)) stop("Alguna de 'variables' no está en 'endog'.")
  nfil <- ceiling(length(sel) / 2); ncol_g <- min(2, length(sel))
  dibujar <- function() {
    op <- par(mfrow = c(nfil, ncol_g), mar = c(4, 4.5, 2.5, 1))
    on.exit(par(op))
    for (j in sel) {
      rango <- range(sapply(regimenes, function(r) resultados[[r]][, j]))
      plot(0:N, rep(0, N + 1), type = "n", ylim = rango,
           xlab = "Horizonte (trimestres)",
           ylab = if (acumular) "Resp. acumulada por 1% de dep." else "Respuesta por 1% de dep.",
           main = endog[j])
      abline(h = 0, col = "grey60")
      for (k in seq_along(regimenes)) lines(0:N, resultados[[regimenes[k]]][, j], col = cols[k], lwd = 2)
      if (j == sel[1]) legend("topleft", bty = "n", lwd = 2, col = cols[seq_along(regimenes)],
                         legend = paste("Régimen", regimenes))
    }
  }

  dibujar()                                            # en el notebook
  .abrir_pdf(paste0(archivo, ".pdf"), width = 5 * ncol_g, height = 4 * nfil)
  dibujar()                                            # en el PDF
  dev.off()

  resumen <- list()
  for (r in regimenes) {
    tab <- data.frame(horizonte = 0:N, round(resultados[[r]][, sel, drop = FALSE], 4),
                      check.names = FALSE)
    cat("\n+++ Régimen", r, "- IRF lineal", if (acumular) "acumulada" else "",
        "por 1% de depreciación inicial (unidades originales) +++\n")
    print(tab, row.names = FALSE)
    resumen[[r]] <- tab
  }
  invisible(list(respuestas = resultados, resumen = resumen))
}


# =====================================================================
# LECTURA "CORTE POR CORTE" DEL ERPT DEPENDIENTE DEL ESTADO
# Cada trimestre del shock tiene su propia trayectoria del ERPT de h = 0 a h = N.
# Estas funciones muestran esas trayectorias en 2D, que es más fácil de leer
# que la superficie 3D. Funcionan con cualquier objeto out_... o inc_...
# (no usan 'res').
# =====================================================================

# --- Referencias: promedio de los trimestres dominados por cada régimen ---
.refs_regimen <- function(out) {
  out <- .normalizar_out(out)
  M <- out$M
  tab <- .ref_por_regimen(out)
  valores <- lapply(seq_len(M), function(m) tab[[m + 1]])
  E <- .nombre_est(out, mayus = TRUE)
  nombres <- if (M == 2) paste(E, c("régimen bajo", "régimen alto"))
             else paste(E, "régimen", seq_len(M))
  list(valores = valores, nombres = nombres, colores = .colores_reg(M))
}

# 1) Trayectorias de trimestres elegidos, superpuestas (gráfico + tabla)
graficar_trayectorias <- function(out, trimestres, archivo = "erpt_trayectorias",
                                  eje_y = "ERPT", titulo = "Trayectoria del ERPT según el trimestre del shock") {
  out <- .normalizar_out(out)
  etq <- .etiqueta_trim(out$anio_dec)
  idx <- match(trimestres, etq)
  if (anyNA(idx)) stop("No encontré: ", paste(trimestres[is.na(idx)], collapse = ", "),
                       ". Usá el formato '1q2002'.")
  N <- out$N
  ref <- .refs_regimen(out)
  nr  <- length(ref$valores)
  cols <- c("black", "darkorange3", "purple", "darkgreen", "brown", "deeppink")[seq_along(idx)]

  dibujar <- function() {
    par(mar = c(4.5, 4.5, 3, 1))
    rango <- range(c(out$erpt[, idx], unlist(ref$valores)), na.rm = TRUE)
    plot(0:N, rep(0, N + 1), type = "n", ylim = rango, xaxt = "n",
         xlab = "Horizonte (trimestres después del shock)", ylab = eje_y, main = titulo)
    axis(1, at = 0:N)
    abline(h = 0, col = "grey70")
    for (m in seq_len(nr)) lines(0:N, ref$valores[[m]], col = ref$colores[m], lwd = 2, lty = 2)
    for (k in seq_along(idx)) lines(0:N, out$erpt[, idx[k]], col = cols[k], lwd = 2.5, type = "b", pch = 19)
    estado <- sapply(idx, function(k) .texto_estado(out, k))
    legend("topleft", bty = "n", cex = 0.85,
           legend = c(paste0(trimestres, " (", estado, ")"), ref$nombres),
           col = c(cols, ref$colores), lwd = c(rep(2.5, length(idx)), rep(2, nr)),
           lty = c(rep(1, length(idx)), rep(2, nr)), pch = c(rep(19, length(idx)), rep(NA, nr)))
  }
  dibujar()
  .abrir_pdf(paste0(archivo, ".pdf"), width = 9, height = 6)
  dibujar()
  dev.off()

  tab <- data.frame(horizonte = 0:N, round(out$erpt[, idx, drop = FALSE], 4), check.names = FALSE)
  colnames(tab)[-1] <- trimestres
  print(tab, row.names = FALSE)
  invisible(tab)
}

# 2) Animación interactiva con barra deslizante (HTML): se mueve el trimestre del shock
graficar_slider <- function(out, archivo = "erpt_slider", eje_y = "ERPT",
                            ms_por_cuadro = 150) {  # milisegundos por trimestre al apretar Play
  out <- .normalizar_out(out)
  N   <- out$N
  M   <- out$M
  K   <- ncol(out$erpt)
  etq <- .etiqueta_trim(out$anio_dec)
  estado <- sapply(seq_len(K), function(k) .texto_estado(out, k))
  df <- data.frame(h = rep(0:N, K),
                   erpt = as.vector(out$erpt),
                   trimestre = rep(etq, each = N + 1),
                   # plotly ordena los cuadros como texto: con un número de 3 dígitos
                   # ("001", "002", ...) el orden alfabético coincide con el cronológico
                   cuadro = rep(sprintf("%03d", seq_len(K)), each = N + 1),
                   color  = rep(.valor_color(out), each = N + 1),
                   estado = rep(estado, each = N + 1))
  df$texto <- paste0(df$trimestre, "<br>h = ", df$h, "<br>", eje_y, " = ",
                     sprintf("%.3f", df$erpt), "<br>", df$estado)
  ref <- .refs_regimen(out)

  g <- plot_ly()
  for (m in seq_along(ref$valores)) {
    g <- g %>% add_lines(data = data.frame(h = 0:N, valor = ref$valores[[m]]), x = ~h, y = ~valor,
                         name = ref$nombres[m], line = list(color = ref$colores[m], dash = "dash"))
  }
  g <- g %>%
    add_trace(data = df, x = ~h, y = ~erpt, frame = ~cuadro, type = "scatter",
              mode = "lines+markers", name = "Trimestre elegido",
              text = ~texto, hoverinfo = "text",
              line = list(color = "grey30", width = 3),
              marker = list(size = 9, color = ~color, cmin = 0, cmax = 1,
                            colorscale = .escala_plotly(M), colorbar = .barra_plotly(M))) %>%
    layout(title = "ERPT según el trimestre del shock (mové la barra o apretá Play)",
           xaxis = list(title = "Horizonte (trimestres después del shock)", dtick = 1),
           yaxis = list(title = eje_y, range = range(c(out$erpt, unlist(ref$valores)), na.rm = TRUE))) %>%
    animation_opts(frame = ms_por_cuadro, transition = 0, redraw = FALSE) %>%
    animation_slider(currentvalue = list(prefix = "Trimestre del shock: "))

  # Reemplazar los números de cuadro por las etiquetas "1q2002" en la barra
  g <- plotly_build(g)
  if (!is.null(g$x$layout$sliders)) {
    pasos <- g$x$layout$sliders[[1]]$steps
    for (i in seq_along(pasos)) {
      n <- suppressWarnings(as.integer(pasos[[i]]$label))
      if (!is.na(n)) g$x$layout$sliders[[1]]$steps[[i]]$label <- etq[n]
    }
  }

  tryCatch(saveWidget(g, paste0(archivo, ".html"), selfcontained = FALSE),
           error = function(e) message("No se pudo guardar el HTML: ", conditionMessage(e)))
  g
}

# 3) GIF animado (para presentaciones): una imagen por trimestre
#    Requiere el paquete gifski: install.packages("gifski")
guardar_gif <- function(out, archivo = "erpt_animacion", cada = 1, segundos_por_cuadro = 0.15,
                        eje_y = "ERPT") {
  if (!requireNamespace("gifski", quietly = TRUE)) {
    stop("Falta el paquete gifski. Instalalo una vez con: install.packages('gifski')")
  }
  out <- .normalizar_out(out)
  N   <- out$N
  etq <- .etiqueta_trim(out$anio_dec)
  ref <- .refs_regimen(out)
  rango <- range(c(out$erpt, unlist(ref$valores)), na.rm = TRUE)
  idx <- seq(1, ncol(out$erpt), by = cada)
  gifski::save_gif({
    for (k in idx) {
      par(mar = c(4.5, 4.5, 3, 1))
      plot(0:N, out$erpt[, k], type = "b", pch = 19, lwd = 3, col = .color_historia(out, k),
           ylim = rango, xaxt = "n", xlab = "Horizonte (trimestres después del shock)", ylab = eje_y,
           main = paste0("Shock en ", etq[k], "   (", .texto_estado(out, k), ")"))
      axis(1, at = 0:N)
      abline(h = 0, col = "grey70")
      for (m in seq_along(ref$valores)) lines(0:N, ref$valores[[m]], col = ref$colores[m], lty = 2, lwd = 2)
      legend("topleft", bty = "n", lty = 2, lwd = 2, col = ref$colores, legend = ref$nombres)
    }
  }, gif_file = paste0(archivo, ".gif"), width = 900, height = 550, delay = segundos_por_cuadro)
  message("GIF guardado en: ", normalizePath(paste0(archivo, ".gif")))
}


# =====================================================================
# ÓPTIMOS LOCALES ALTERNATIVOS
# Requiere que 'res' contenga el modelo sin identificar (res$tvar_model):
# ver el cambio indicado en sstvars_wrapper.r.
# =====================================================================

# 1) Tabla de las mejores soluciones encontradas en las rondas de estimación
#    (ordenadas de mayor a menor log-verosimilitud, sin repetidas)
explorar_optimos <- function(res, thr_var = "ipc", n = 30, regimen_alto = 2) {
  mod <- res$tvar_model
  if (is.null(mod)) stop("'res' no contiene tvar_model. Actualizá sstvars_wrapper.r ",
                         "y volvé a correr run_vlstar() con load_presaved = TRUE.")
  sd_z <- unname(res$scaling_sd[thr_var])
  n <- min(n, length(mod$all_logliks))

  filas <- lapply(seq_len(n), function(j) {
    avisos <- character(0)
    alt <- withCallingHandlers(
      tryCatch(alt_stvar(mod, which_largest = j, calc_std_errors = FALSE),
               error = function(e) NULL),
      warning = function(w) {
        avisos <<- c(avisos, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    if (is.null(alt)) return(NULL)
    ll <- tryCatch(as.numeric(if (!is.null(alt$loglik)) alt$loglik else logLik(alt)),
                   error = function(e) NA)
    tipo_j <- .tipo_transicion(alt)
    wp <- .param_transicion(alt)
    a  <- alt$transition_weights[, min(regimen_alto, ncol(alt$transition_weights))]
    dom <- apply(alt$transition_weights, 1, which.max)
    es_log <- identical(tipo_j, "logistic")
    data.frame(which_largest    = j,
               loglik           = round(ll, 2),
               umbral_trim_pct  = if (es_log) as.character(round(100 * (exp(wp[1] * sd_z) - 1), 1))
                                  else if (identical(tipo_j, "threshold"))
                                    paste(round(100 * (exp(sort(wp) * sd_z) - 1), 1), collapse = " / ")
                                  else NA,
               gamma            = if (es_log) round(wp[2], 2) else NA,
               gl_student       = round(.grados_libertad(alt), 1),
               trim_por_regimen = paste(tabulate(dom, alt$model$M), collapse = "/"),
               trim_reg_alto    = sum(a > 0.9),
               trim_transicion  = sum(a > 0.1 & a < 0.9),
               raiz_casi_unit   = any(grepl("unit", avisos)),
               cov_casi_singular = any(grepl("singular", avisos)),
               descartar        = length(avisos) > 0)
  })
  tab <- do.call(rbind, filas)
  tab <- tab[!duplicated(tab[, c("loglik", "umbral_trim_pct", "gamma")]), ]
  cat("umbral_trim_pct: umbral(es) en % trimestral (logístico y TVAR; vacío en mlogit) | gamma: solo logístico | gl_student: grados de libertad (solo errores t)\n",
      "trim_por_regimen: trimestres en que domina cada régimen (1/2/...)\n",
      "raiz_casi_unit / cov_casi_singular / descartar: TRUE = sstvars avisó un problema; no usar ese óptimo\n\n")
  print(tab, row.names = FALSE)
  invisible(tab)
}

# 2) Construir un 'res' alternativo con la solución elegida, identificado por
#    heteroscedasticidad como en el wrapper. Funciona con erpt_historia(),
#    graficar_girf_original(), graficar_irf_original(), graficar_transicion()
#    y la celda de gamma (no con print_passthrough()).
usar_optimo <- function(res, which_largest) {
  alt <- alt_stvar(res$tvar_model, which_largest = which_largest, calc_std_errors = FALSE)
  ident <- if (.dist_errores(alt) %in% c("ind_Student", "ind_skewed_t")) "non-Gaussianity" else "heteroskedasticity"
  est <- fitSSTVAR(stvar = alt, identification = ident)
  list(fit_struct = est, tvar_model = alt, scaling_sd = res$scaling_sd,
       which_largest = which_largest)
}


# =====================================================================
# CRITERIOS DE INFORMACIÓN para comparar modelos con distinto número de regímenes
# Usa el modelo SIN identificar (forma reducida), para que la comparación
# no dependa de las restricciones de la identificación por heteroscedasticidad.
# Menor valor = mejor modelo. Requiere res$tvar_model (wrapper v2).
# =====================================================================
criterios_info <- function(res, etiqueta) {
  m <- res$tvar_model
  if (is.null(m)) stop("'res' no contiene tvar_model. Cargá sstvars_wrapper_v2.r ",
                       "y volvé a correr el modelo con load_presaved = TRUE.")
  ll <- tryCatch(as.numeric(if (!is.null(m$loglik)) m$loglik else logLik(m)),
                 error = function(e) NA)
  k  <- length(m$params)
  Tn <- nrow(m$data) - m$model$p
  data.frame(modelo     = etiqueta,
             regimenes  = m$model$M,
             parametros = k,
             loglik     = round(ll, 2),
             AIC        = round(-2 * ll + 2 * k, 1),
             HQ         = round(-2 * ll + 2 * k * log(log(Tn)), 1),
             BIC        = round(-2 * ll + k * log(Tn), 1),
             check.names = FALSE)
}


# =====================================================================
# LECTURA DEL MODELO EN UNIDADES ORIGINALES (2 o 3 regímenes)
# Umbrales en % de inflación trimestral y anualizada, trimestres por régimen
# y promedios observados de cada variable en los trimestres de cada régimen.
# =====================================================================
resumen_modelo <- function(res, data, endog, trans_var = "ipc", modelo = NULL, escalas = NULL) {
  if (is.null(modelo)) {
    modelo <- .buscar(res, function(x) inherits(x, "stvar") &&
                        !identical(x$model$identification, "reduced_form"))
  }
  if (is.null(modelo)) stop("No encontré el modelo estructural en 'res'. Pasalo con modelo = ...")
  if (is.null(escalas)) escalas <- .buscar_nombre(res, "scaling_sd")
  if (is.null(escalas)) escalas <- sapply(data[, endog], sd, na.rm = TRUE)

  M    <- modelo$model$M
  p    <- modelo$model$p
  tipo <- .tipo_transicion(modelo)
  sd_z <- unname(escalas[match(trans_var, endog)])
  pct_t <- function(x) round(100 * (exp(x) - 1), 1)
  pct_a <- function(x) round(100 * (exp(4 * x) - 1), 0)

  cat("Modelo:", tipo, "| regímenes:", M, "| rezagos:", p,
      "| errores:", .dist_errores(modelo), "\n")
  if (!is.na(.grados_libertad(modelo)))
    cat("Grados de libertad de la t de Student:", round(.grados_libertad(modelo), 2), "\n")
  if (identical(tipo, "threshold")) {
    u <- sort(.param_transicion(modelo)) * sd_z
    cat("Umbrales:", paste0(pct_t(u), "% trimestral (", pct_a(u), "% anualizada)",
                             collapse = " | "), "\n")
  } else if (identical(tipo, "mlogit")) {
    idx_z <- match(trans_var, endog)
    tz    <- (p + 1):nrow(modelo$data)
    z     <- as.matrix(modelo$data)[tz - 1, idx_z] * sd_z
    cr    <- .cruces_empiricos(z, modelo$transition_weights[tz - p, , drop = FALSE])
    cat("Transición mlogit (parámetros gamma, escalados):", paste(round(.param_transicion(modelo), 3), collapse = ", "), "\n")
    cat("Cambios de régimen dominante según la inflación de t-1:\n")
    if (nrow(cr)) print(data.frame(de_regimen = cr$de, a_regimen = cr$a,
                                   inflacion_trim_pct = pct_t(cr$z),
                                   inflacion_anual_pct = pct_a(cr$z)), row.names = FALSE)
    if (nrow(cr) > M - 1) cat("ATENCIÓN: más cambios que regímenes - 1: los regímenes no están ordenados por inflación.\n")
  } else {
    wp <- .param_transicion(modelo)
    cat("Umbral (c):", pct_t(wp[1] * sd_z), "% trimestral (", pct_a(wp[1] * sd_z),
        "% anualizada) | gamma:", round(wp[2], 2), "\n")
  }

  T_obs <- nrow(modelo$data)
  t_idx <- (p + 1):T_obs
  pesos <- modelo$transition_weights[t_idx - p, , drop = FALSE]
  reg   <- apply(pesos, 1, which.max)
  cat("Trimestres por régimen (régimen de mayor peso):",
      paste0("régimen ", seq_len(M), " = ", tabulate(reg, M), collapse = " | "), "\n\n")

  prom <- sapply(seq_len(M), function(m) colMeans(data[t_idx[reg == m], endog, drop = FALSE], na.rm = TRUE))
  prom <- round(100 * (exp(prom) - 1), 2)
  colnames(prom) <- paste("Régimen", seq_len(M))
  cat("Promedio observado de cada variable en los trimestres de cada régimen\n",
      "(variación trimestral en %, a partir de las diferencias logarítmicas):\n", sep = "")
  print(prom)
  invisible(list(regimen = reg, promedios = prom))
}

# =====================================================================
# IDENTIFICAR EL SHOCK CAMBIARIO
# Imprime el efecto de impacto de cada shock (unidades originales) y sugiere
# como cambiario al que concentra su efecto sobre E. Confirmar mirando la tabla:
# el cambiario mueve mucho a E y muy poco al ipc en el impacto.
# =====================================================================
identificar_shock_cambiario <- function(res, endog, shock_var = "E", resp_var = "ipc",
                                        modelo = NULL, escalas = NULL) {
  if (is.null(modelo)) {
    modelo <- .buscar(res, function(x) inherits(x, "stvar") &&
                        !identical(x$model$identification, "reduced_form"))
  }
  if (is.null(modelo)) stop("No encontré el modelo estructural en 'res'. Pasalo con modelo = ...")
  if (is.null(escalas)) escalas <- .buscar_nombre(res, "scaling_sd")
  if (is.null(escalas)) stop("No encontré los desvíos de escalado en 'res'. Pasalos con escalas = ...")
  d <- length(endog)
  M <- modelo$model$M
  ie <- match(shock_var, endog); ip <- match(resp_var, endog)
  no_gauss <- .dist_errores(modelo) %in% c("ind_Student", "ind_skewed_t")

  # Efecto de impacto en cada régimen. Con identificación por heteroscedasticidad, el
  # impacto de un mismo shock es proporcional entre regímenes; con no normalidad
  # (ind_Student) cada régimen tiene su propia matriz de impacto B_m.
  imp_reg <- lapply(seq_len(M), function(m) linear_IRF(modelo, N = 1, regime = m)$point_est[, , 1])
  part_reg <- sapply(imp_reg, function(x) abs(x[ie, ]) / sqrt(colSums(x^2)))   # [shock, régimen]
  part_reg <- matrix(part_reg, nrow = d)
  sug_reg  <- apply(part_reg, 2, which.max)

  imp_esc <- imp_reg[[1]]
  imp <- imp_esc * unname(escalas)                                        # unidades originales
  dimnames(imp) <- list(respuesta = endog, shock = paste0("shock_", seq_len(d)))
  part_e <- part_reg[, 1]
  sugerido <- which.max(part_e)

  cat("Efecto de impacto de cada shock", if (M > 1 && no_gauss) "(régimen 1)" else "", "(unidades originales):\n")
  print(round(imp, 4))
  cat("\nParticipación de", shock_var, "en el efecto de cada shock",
      if (M > 1 && no_gauss) "(régimen 1)" else "", ":\n")
  print(round(setNames(part_e, colnames(imp)), 2))
  cat("\nERPT de impacto de cada shock (", resp_var, "/", shock_var, "):\n", sep = "")
  print(round(imp[ip, ] / imp[ie, ], 3))

  if (M > 1 && no_gauss) {
    for (m in 2:M) {
      im <- imp_reg[[m]] * unname(escalas)
      dimnames(im) <- dimnames(imp)
      cat("\n--- Régimen", m, "---\nEfecto de impacto (unidades originales):\n"); print(round(im, 4))
      cat("Participación de", shock_var, ":\n"); print(round(setNames(part_reg[, m], colnames(imp)), 2))
      cat("ERPT de impacto:\n"); print(round(im[ip, ] / im[ie, ], 3))
    }
    signos <- sapply(imp_reg, function(x) sign(x[ie, sug_reg[1]]))
    if (length(unique(signos)) > 1) {
      cat("\n!!! ATENCIÓN: el shock", sug_reg[1], "mueve", shock_var, "con signo distinto según el régimen:",
          paste0("régimen ", seq_len(M), " -> ", ifelse(signos > 0, "+", "-"), collapse = " | "),
          "\n    El signo de las columnas de B_2, ..., B_M está débilmente identificado: un mismo shock no debería",
          "\n    depreciar en un régimen y apreciar en otro. No calcular el ERPT hasta estimar el modelo con el",
          "\n    signo corregido (fitSSTVAR con B_pm_reg, B_perm y B_signs) y comparar la verosimilitud.\n")
    }
    if (length(unique(sug_reg)) > 1) {
      cat("\n!!! ATENCIÓN: el shock con mayor participación de", shock_var, "no es el mismo en todos los regímenes:",
          paste0("régimen ", seq_len(M), " -> shock ", sug_reg, collapse = " | "),
          "\n    Con identificación por no normalidad el orden de las columnas de B puede diferir entre regímenes.",
          "\n    No calcular el ERPT hasta revisar el orden de los shocks (ver reorder_B_columns / fitSSTVAR con B_perm).\n")
    } else {
      cat("\nEl mismo shock concentra el efecto sobre", shock_var, "en todos los regímenes.\n")
    }
  }
  cat("\n>>> Shock sugerido como cambiario:", sugerido,
      "(confirmalo: debe mover mucho a", shock_var, "y poco a", resp_var, ")\n")
  invisible(sugerido)
}


# =====================================================================
# TABLA DE MODELOS GUARDADOS con criterios de información
# Recorre la carpeta models/ SIN reestimar ni recalcular GIRF, y lista los
# modelos (forma reducida) estimados con exactamente los datos de 'data'
# (por ejemplo, df_ERPT_Arg para "sin cappeo" o df_ERPT_Arg_capped para "con cappeo").
# Menor AIC / HQ / BIC = mejor modelo.
# Ojo: T_efectivo cambia con los rezagos (se pierden p observaciones), así que la
# comparación entre distintos rezagos es aproximada.
# =====================================================================
tabla_modelos_guardados <- function(data, carpeta = "models") {
  archivos <- list.files(carpeta, pattern = "\\.rds$", full.names = TRUE)
  archivos <- archivos[!grepl("^structured_", basename(archivos))]
  filas <- lapply(archivos, function(f) {
    m <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(m) || is.null(m$data) || is.null(m$model)) return(NULL)
    Y <- as.matrix(m$data)
    vars <- colnames(Y)
    if (is.null(vars) || !all(vars %in% names(data))) return(NULL)
    X <- as.matrix(data[, vars])
    if (nrow(X) != nrow(Y)) return(NULL)
    X <- sweep(X, 2, apply(X, 2, sd, na.rm = TRUE), "/")
    if (!isTRUE(all.equal(unname(X), unname(Y), tolerance = 1e-6))) return(NULL)
    ll <- tryCatch(as.numeric(if (!is.null(m$loglik)) m$loglik else logLik(m)),
                   error = function(e) NA)
    k  <- length(m$params)
    p  <- m$model$p
    Tn <- nrow(Y) - p
    data.frame(commodity  = vars[1],
               modelo     = if (grepl("^vlstar", basename(f))) "VLSTAR" else "TVAR",
               regimenes  = m$model$M,
               rezagos    = p,
               parametros = k,
               T_efectivo = Tn,
               loglik     = round(ll, 2),
               AIC        = round(-2 * ll + 2 * k, 1),
               HQ         = round(-2 * ll + 2 * k * log(log(Tn)), 1),
               BIC        = round(-2 * ll + k * log(Tn), 1),
               check.names = FALSE)
  })
  tab <- do.call(rbind, filas)
  if (is.null(tab)) { cat("No encontré modelos estimados con esos datos.\n"); return(invisible(NULL)) }
  tab <- tab[order(tab$commodity, tab$modelo, tab$regimenes, tab$rezagos), ]
  print(tab, row.names = FALSE)
  invisible(tab)
}


# =====================================================================
# 3D: FECHA DEL SHOCK x RÉGIMEN x ERPT a un horizonte fijo (por defecto h = 8)
# - Con 2 regímenes (VLSTAR), el eje del régimen es alfa (peso del régimen alto, 0 a 1).
# - Con 3 regímenes (TVAR), el eje del régimen es 1, 2 o 3.
# Usar con out_..._acum para el ERPT relativo a la depreciación acumulada.
# No usa 'res'.
# =====================================================================
graficar_erpt_regimen_3d <- function(out, h = 8, archivo = "erpt_regimen_3d",
                                     eje_z = "ERPT", unir = TRUE,
                                     titulo = NULL) {
  out <- .normalizar_out(out)
  M <- out$M
  if (h > out$N) stop("El horizonte h no puede ser mayor que ", out$N, ".")
  x  <- out$anio_dec
  yv <- if (M == 2) out$alpha else out$regimen
  z  <- out$erpt[h + 1, ]
  etq <- .etiqueta_trim(x)
  eje_y <- if (M == 2) "\u03b1 régimen alto" else "Régimen"
  if (is.null(titulo)) titulo <- paste0(eje_z, " a h = ", h, " según fecha del shock y régimen")
  tk <- .ticks_trim(x)
  vc <- .valor_color(out)

  # ---- Tabla resumen por régimen ----
  if (M == 2) {
    grupos <- list("Régimen bajo (alfa < 0,1)" = out$alpha < 0.1,
                   "Transición" = out$alpha >= 0.1 & out$alpha <= 0.9,
                   "Régimen alto (alfa > 0,9)" = out$alpha > 0.9)
  } else {
    grupos <- setNames(lapply(seq_len(M), function(m) out$regimen == m),
                       paste("Régimen", seq_len(M)))
  }
  resumen <- data.frame(grupo = names(grupos),
                        trimestres = sapply(grupos, sum),
                        promedio = round(sapply(grupos, function(g) mean(z[g])), 3),
                        mediana  = round(sapply(grupos, function(g) if (any(g)) median(z[g]) else NA), 3),
                        minimo   = round(sapply(grupos, function(g) if (any(g)) min(z[g]) else NA), 3),
                        maximo   = round(sapply(grupos, function(g) if (any(g)) max(z[g]) else NA), 3),
                        row.names = NULL)
  cat(eje_z, "a h =", h, "por régimen:\n")
  print(resumen, row.names = FALSE)

  # ---- 3D interactivo ----
  texto <- paste0(etq, "<br>", .letra_alfa_html(M, out, seq_along(x)),
                  "<br>", eje_z, " (h = ", h, "): ", sprintf("%.3f", z))
  g <- plot_ly(x = x, y = yv, z = z, type = "scatter3d",
               mode = if (unir) "lines+markers" else "markers",
               text = texto, hoverinfo = "text",
               line = list(color = "rgba(120,120,120,0.5)", width = 2),
               marker = list(size = 4, color = vc, cmin = 0, cmax = 1,
                             colorscale = .escala_plotly(M), colorbar = .barra_plotly(M))) %>%
    layout(title = titulo,
           scene = list(xaxis = list(title = "Trimestre del shock", tickvals = tk$at, ticktext = tk$lab),
                        yaxis = if (M == 2) list(title = eje_y, range = c(0, 1))
                                else list(title = eje_y, tickvals = seq_len(M), range = c(0.7, M + 0.3)),
                        zaxis = list(title = paste0(eje_z, " (h = ", h, ")"))))

  # ---- 3D estático para la tesis ----
  ylim <- if (M == 2) c(0, 1) else c(0.7, M + 0.3)
  zlim <- range(z, na.rm = TRUE)
  cols <- sapply(seq_along(x), function(k) .color_historia(out, k))
  dibujar <- function() {
    pm <- persp(x = range(x), y = ylim, z = matrix(zlim[1], 2, 2), zlim = zlim,
                theta = -35, phi = 22, expand = 0.7, ticktype = "detailed",
                col = NA, border = NA, box = TRUE,
                xlab = "Trimestre del shock", ylab = if (M == 2) "alfa régimen alto" else "Régimen",
                zlab = paste0(eje_z, " (h = ", h, ")"), main = titulo, cex.lab = 0.8, cex.axis = 0.7)
    base <- trans3d(x, yv, rep(zlim[1], length(x)), pm)
    pts  <- trans3d(x, yv, z, pm)
    segments(base$x, base$y, pts$x, pts$y, col = adjustcolor(cols, 0.25))
    if (unir) lines(pts, col = "grey60")
    points(pts, pch = 19, cex = 0.6, col = cols)
  }
  .abrir_pdf(paste0(archivo, ".pdf"), width = 9, height = 7)
  dibujar()
  dev.off()

  tryCatch(saveWidget(g, paste0(archivo, ".html"), selfcontained = FALSE),
           error = function(e) message("No se pudo guardar el HTML: ", conditionMessage(e)))
  g
}

# Texto del estado para el cuadro del mouse
.letra_alfa_html <- function(M, out, k) {
  if (M == 2) paste0("\u03b1 = ", sprintf("%.2f", out$alpha[k]))
  else paste0("Régimen ", out$regimen[k])
}


# =====================================================================
# RED 3D CONTINUA: FECHA DEL SHOCK x RÉGIMEN x ERPT a un horizonte fijo
# Cada régimen es una "cinta" (superficie) a lo largo del tiempo, cuya altura es
# el ERPT a h trimestres. La cinta existe solo en los trimestres en que la
# economía estuvo en ese régimen: los huecos muestran los cambios de régimen.
# - 3 regímenes (TVAR): cintas en régimen 1, 2 y 3.
# - 2 regímenes (VLSTAR): cintas en régimen bajo (alfa < 0,1), transición y alto (alfa > 0,9).
# No usa 'res'.
# =====================================================================
graficar_erpt_regimen_red <- function(out, h = 8, archivo = "erpt_regimen_red",
                                      eje_z = "ERPT", titulo = NULL,
                                      ancho_cinta = 0.35, cada_anios_red = 2) {
  out <- .normalizar_out(out)
  M <- out$M
  if (h > out$N) stop("El horizonte h no puede ser mayor que ", out$N, ".")
  x   <- out$anio_dec
  z   <- out$erpt[h + 1, ]
  etq <- .etiqueta_trim(x)
  K   <- length(x)

  # Grupos (filas del eje "régimen")
  if (M == 2) {
    grupo <- ifelse(out$alpha < 0.1, 1, ifelse(out$alpha > 0.9, 3, 2))
    nombres <- c("Bajo", "Transición", "Alto")
    cols <- c("steelblue", "orange", "firebrick")
  } else {
    grupo <- out$regimen
    nombres <- paste("Régimen", seq_len(M))
    cols <- .colores_reg(M)
  }
  G <- length(nombres)
  if (is.null(titulo)) titulo <- paste0(eje_z, " a h = ", h, " según fecha del shock y régimen")

  # Cada trimestre ocupa una franja de un trimestre de ancho (dos columnas de la red),
  # así los trimestres consecutivos en el mismo régimen quedan unidos
  # y un trimestre aislado igual se ve.
  xx <- as.vector(rbind(x - 0.12, x + 0.12))
  yy <- as.vector(rbind(seq_len(G) - ancho_cinta, seq_len(G) + ancho_cinta))
  Z  <- matrix(NA_real_, nrow = length(yy), ncol = length(xx))
  C  <- matrix(NA_real_, nrow = length(yy), ncol = length(xx))
  TX <- matrix("", nrow = length(yy), ncol = length(xx))
  for (k in seq_len(K)) {
    g  <- grupo[k]
    fi <- c(2 * g - 1, 2 * g)
    co <- c(2 * k - 1, 2 * k)
    Z[fi, co]  <- z[k]
    C[fi, co]  <- (g - 1) / (G - 1)
    TX[fi, co] <- paste0(etq[k], "<br>", nombres[g], "<br>", eje_z, " (h = ", h, "): ",
                         sprintf("%.3f", z[k]))
  }

  # Resumen por régimen
  resumen <- data.frame(grupo = nombres,
                        trimestres = tabulate(grupo, G),
                        promedio = round(sapply(seq_len(G), function(g) mean(z[grupo == g])), 3),
                        mediana  = round(sapply(seq_len(G), function(g) if (any(grupo == g)) median(z[grupo == g]) else NA), 3),
                        minimo   = round(sapply(seq_len(G), function(g) if (any(grupo == g)) min(z[grupo == g]) else NA), 3),
                        maximo   = round(sapply(seq_len(G), function(g) if (any(grupo == g)) max(z[grupo == g]) else NA), 3))
  cat(eje_z, "a h =", h, "por régimen:\n")
  print(resumen, row.names = FALSE)

  # Escala de colores discreta por régimen
  esc <- list()
  for (g in seq_len(G)) {
    esc[[length(esc) + 1]] <- list((g - 1) / G, cols[g])
    esc[[length(esc) + 1]] <- list(g / G, cols[g])
  }
  tk <- .ticks_trim(x)
  lim_x <- range(xx)
  g3d <- plot_ly(x = xx, y = yy, z = Z, type = "surface",
                 surfacecolor = C, cmin = 0, cmax = 1, colorscale = esc,
                 text = TX, hoverinfo = "text", connectgaps = FALSE,
                 colorbar = list(title = "Régimen", tickvals = (seq_len(G) - 0.5) / G,
                                 ticktext = nombres),
                 contours = list(
                   x = list(show = TRUE, start = floor(lim_x[1]), end = ceiling(lim_x[2]),
                            size = cada_anios_red, color = "rgba(40,40,40,0.5)", width = 1),
                   y = list(show = TRUE, start = 1 - ancho_cinta, end = G + ancho_cinta,
                            size = ancho_cinta, color = "rgba(40,40,40,0.5)", width = 1))) %>%
    layout(title = titulo,
           scene = list(xaxis = list(title = "Trimestre del shock", tickvals = tk$at, ticktext = tk$lab),
                        yaxis = list(title = "Régimen", tickvals = seq_len(G), ticktext = nombres),
                        zaxis = list(title = paste0(eje_z, " (h = ", h, ")")),
                        aspectratio = list(x = 2, y = 1, z = 0.8)))

  # Versión estática (PDF) con persp: misma red, coloreada por régimen
  .abrir_pdf(paste0(archivo, ".pdf"), width = 11, height = 7)
  nf <- length(yy) - 1; nc <- length(xx) - 1
  col_fac <- matrix(NA_character_, nf, nc)
  for (i in seq_len(nf)) for (j in seq_len(nc)) {
    v <- C[i:(i + 1), j:(j + 1)]
    if (all(!is.na(v))) col_fac[i, j] <- cols[round(v[1] * (G - 1)) + 1]
  }
  persp(x = xx, y = yy, z = t(Z), theta = -35, phi = 25, expand = 0.6,
        col = t(col_fac), border = NA, shade = 0.3, ticktype = "detailed",
        xlab = "Trimestre del shock", ylab = "Régimen", zlab = paste0(eje_z, " (h = ", h, ")"),
        main = titulo, cex.lab = 0.8, cex.axis = 0.7)
  dev.off()

  tryCatch(saveWidget(g3d, paste0(archivo, ".html"), selfcontained = FALSE),
           error = function(e) message("No se pudo guardar el HTML: ", conditionMessage(e)))
  g3d
}


# =====================================================================
# CORTINA 3D CONTINUA (sin huecos, solo datos observados)
# Una única superficie a lo largo del tiempo:
#   - eje X: fecha del shock
#   - eje Y: estado del período (alfa con 2 regímenes; régimen 1-2-3 con 3 regímenes)
#   - eje Z: ERPT a h trimestres (el borde superior de la cortina)
# Cuando la economía cambia de régimen, la cortina se desplaza en el eje Y.
# No usa 'res'.
# =====================================================================
graficar_erpt_cortina <- function(out, h = 8, archivo = "erpt_cortina",
                                  eje_z = "ERPT", titulo = NULL, cada_anios_red = 2,
                                  theta = -35, phi = 25) {
  out <- .normalizar_out(out)
  M <- out$M
  if (h > out$N) stop("El horizonte h no puede ser mayor que ", out$N, ".")
  x   <- out$anio_dec
  yv  <- if (M == 2) out$alpha else out$regimen
  z   <- out$erpt[h + 1, ]
  K   <- length(x)
  etq <- .etiqueta_trim(x)
  vc  <- .valor_color(out)
  base <- min(0, min(z, na.rm = TRUE))
  eje_y <- if (M == 2) "\u03b1 régimen alto" else "Régimen"
  if (is.null(titulo)) titulo <- paste0(eje_z, " a h = ", h, " según fecha del shock y régimen")

  estado <- sapply(seq_len(K), function(k) .texto_estado(out, k))
  cat(eje_z, "a h =", h, "-", .nombre_est(out), "por régimen:\n")
  print(round(.ref_por_regimen(out, matrix(z, nrow = 1), 0)[, -1, drop = FALSE], 3),
        row.names = FALSE)

  # Superficie de 2 filas (base y borde superior) x K columnas (trimestres)
  X  <- rbind(x, x)
  Y  <- rbind(yv, yv)
  Z  <- rbind(rep(base, K), z)
  C  <- rbind(vc, vc)
  TX <- rbind(paste0(etq, "<br>", estado, "<br>", eje_z, " (h = ", h, "): ", sprintf("%.3f", z)),
              paste0(etq, "<br>", estado, "<br>", eje_z, " (h = ", h, "): ", sprintf("%.3f", z)))
  tk <- .ticks_trim(x)
  lim_x <- range(x)

  g <- plot_ly() %>%
    add_surface(x = X, y = Y, z = Z, surfacecolor = C, cmin = 0, cmax = 1,
                colorscale = .escala_plotly(M), colorbar = .barra_plotly(M),
                text = TX, hoverinfo = "text", opacity = 0.9,
                contours = list(x = list(show = TRUE, start = floor(lim_x[1]),
                                         end = ceiling(lim_x[2]), size = cada_anios_red,
                                         color = "rgba(40,40,40,0.6)", width = 1),
                                z = list(show = TRUE, color = "rgba(40,40,40,0.3)", width = 1))) %>%
    add_trace(x = x, y = yv, z = z, type = "scatter3d", mode = "lines",
              line = list(color = "black", width = 3), text = TX[2, ], hoverinfo = "text",
              showlegend = FALSE) %>%
    layout(title = titulo,
           scene = list(xaxis = list(title = "Trimestre del shock", tickvals = tk$at, ticktext = tk$lab),
                        yaxis = if (M == 2) list(title = eje_y, range = c(0, 1))
                                else list(title = eje_y, tickvals = seq_len(M), range = c(0.8, M + 0.2)),
                        zaxis = list(title = paste0(eje_z, " (h = ", h, ")")),
                        aspectratio = list(x = 2.2, y = 1, z = 0.8)))

  # Versión estática (PDF)
  cols <- sapply(seq_len(K), function(k) .color_historia(out, k))
  ylim <- if (M == 2) c(0, 1) else c(0.8, M + 0.2)
  .abrir_pdf(paste0(archivo, ".pdf"), width = 11, height = 7)
  pm <- persp(x = range(x), y = ylim, z = matrix(base, 2, 2), zlim = range(c(base, z)),
              theta = theta, phi = phi, expand = 0.6, ticktype = "detailed",
              col = NA, border = NA, box = TRUE,
              xlab = "Trimestre del shock", ylab = if (M == 2) "alfa régimen alto" else "Régimen",
              zlab = paste0(eje_z, " (h = ", h, ")"), main = titulo, cex.lab = 0.8, cex.axis = 0.7)
  for (k in seq_len(K - 1)) {
    p1 <- trans3d(c(x[k], x[k + 1], x[k + 1], x[k]),
                  c(yv[k], yv[k + 1], yv[k + 1], yv[k]),
                  c(base, base, z[k + 1], z[k]), pm)
    polygon(p1, col = adjustcolor(cols[k], 0.85), border = adjustcolor("grey20", 0.3), lwd = 0.3)
  }
  lines(trans3d(x, yv, z, pm), col = "black", lwd = 1)
  dev.off()

  tryCatch(saveWidget(g, paste0(archivo, ".html"), selfcontained = FALSE),
           error = function(e) message("No se pudo guardar el HTML: ", conditionMessage(e)))
  g
}


# =====================================================================
# 2D: ERPT a horizontes fijos a lo largo del tiempo, con BANDAS DE COLOR
# de fondo según el régimen de cada trimestre.
#   - 3 regímenes: azul / naranja / rojo.
#   - 2 regímenes: degradé de azul (alfa = 0) a rojo (alfa = 1).
# No usa 'res'.
# =====================================================================
graficar_erpt_2d_bandas <- function(out, horizontes = 8, archivo = "erpt_2d_bandas",
                                    eje_y = "ERPT", titulo = NULL,
                                    promedios = TRUE, intensidad = 0.30,
                                    ancho = 11, alto = 5.5) {
  out <- .normalizar_out(out)
  M <- out$M
  if (any(horizontes > out$N)) stop("Algún horizonte es mayor que ", out$N, ".")
  x  <- out$anio_dec
  K  <- length(x)
  tk <- .ticks_trim(x)
  if (is.null(titulo)) {
    titulo <- if (length(horizontes) == 1) paste0(eje_y, " a h = ", horizontes, " según el trimestre del shock")
              else paste0(eje_y, " a horizontes seleccionados según el trimestre del shock")
  }
  cols_banda <- sapply(seq_len(K), function(k) adjustcolor(.color_historia(out, k), intensidad))
  cols_lin <- c("black", "grey40", "darkgreen", "purple")[seq_along(horizontes)]
  lty_lin  <- c(1, 2, 3, 4)[seq_along(horizontes)]

  # Promedio (o mediana) del ERPT por régimen, para cada horizonte
  tab <- .ref_por_regimen(out)[horizontes + 1, , drop = FALSE]
  cat(.nombre_est(out, mayus = TRUE), "de", eje_y, "por régimen (trimestres con peso > 0.9):\n")
  print(round(tab, 3), row.names = FALSE)

  dibujar <- function() {
    op <- par(mar = c(4.5, 4.5, 3, 1)); on.exit(par(op))
    rango <- range(out$erpt[horizontes + 1, ], na.rm = TRUE)
    rango <- rango + c(-0.03, 0.03) * diff(rango)
    plot(NA, xlim = range(x) + c(-0.125, 0.125), ylim = rango, xaxt = "n", xaxs = "i",
         xlab = "Trimestre del shock", ylab = eje_y, main = titulo)
    rect(x - 0.125, rango[1] - 1, x + 0.125, rango[2] + 1, col = cols_banda, border = NA)
    axis(1, at = tk$at, labels = tk$lab)
    abline(h = 0, col = "grey60")
    # Promedio por régimen (líneas punteadas del color del régimen), solo con un horizonte
    if (promedios && length(horizontes) == 1) {
      cr <- .colores_reg(M)
      for (m in seq_len(M)) {
        v <- tab[1, m + 1]
        if (!is.na(v)) abline(h = v, col = cr[m], lty = 3, lwd = 1.5)
      }
    }
    for (i in seq_along(horizontes)) {
      lines(x, out$erpt[horizontes[i] + 1, ], col = cols_lin[i], lty = lty_lin[i], lwd = 1.6)
    }
    # Leyendas
    if (M == 2) {
      leg_reg <- c("Régimen bajo (\u03b1 = 0)", "Régimen alto (\u03b1 = 1)")
      if (!isTRUE(capabilities("cairo"))) leg_reg <- c("Régimen bajo (alfa = 0)", "Régimen alto (alfa = 1)")
      col_reg <- .colores_reg(2)
    } else {
      leg_reg <- paste("Régimen", seq_len(M)); col_reg <- .colores_reg(M)
    }
    legend("topleft", bty = "n", cex = 0.85,
           legend = c(paste0("h = ", horizontes, " trimestres"), leg_reg),
           col = c(cols_lin, adjustcolor(col_reg, intensidad + 0.2)),
           lty = c(lty_lin, rep(NA, length(leg_reg))), lwd = c(rep(1.6, length(horizontes)), rep(NA, length(leg_reg))),
           pch = c(rep(NA, length(horizontes)), rep(15, length(leg_reg))), pt.cex = 2)
    if (promedios && length(horizontes) == 1)
      mtext(paste0("Líneas punteadas: ", .nombre_est(out), " de cada régimen"), side = 3, line = 0.2, cex = 0.75, adj = 1)
  }

  dibujar()                                       # en el notebook
  .abrir_pdf(paste0(archivo, ".pdf"), width = ancho, height = alto)
  dibujar()                                       # en el PDF
  dev.off()
  invisible(tab)
}


# =====================================================================
# MANTA DE ALAMBRE: FECHA DEL SHOCK x HORIZONTE x ERPT
# Solo valores observados: cada trimestre con el régimen que realmente existió.
# La red se colorea según ese régimen. 'suavizar' (opcional) aplica un promedio
# móvil centrado de k trimestres sobre la fecha (modifica los valores: aclararlo).
# No usa 'res'.
# =====================================================================
graficar_erpt_manta <- function(out, archivo = "erpt_manta", eje_z = "ERPT",
                                titulo = NULL, suavizar = 1, cada = 1,
                                theta = -40, phi = 28, fondo = "white",
                                ancho = 11, alto = 8, opaca = TRUE) {
  out <- .normalizar_out(out)
  N <- out$N
  Z <- out$erpt                                  # [horizonte, trimestre]
  if (suavizar > 1) {
    Z <- t(apply(Z, 1, function(v) {
      s <- stats::filter(v, rep(1 / suavizar, suavizar), sides = 2)
      ifelse(is.na(s), v, s)
    }))
  }
  idx <- seq(1, ncol(Z), by = cada)
  x   <- out$anio_dec[idx]
  Z   <- Z[, idx, drop = FALSE]
  K   <- length(x)
  if (is.null(titulo)) titulo <- if (suavizar > 1)
    paste0(eje_z, " dependiente del estado (promedio móvil de ", suavizar, " trimestres)")
    else paste0(eje_z, " dependiente del estado")
  col_lin <- sapply(idx, function(k) .color_historia(out, k))

  # ---- PDF / pantalla: red de alambre con persp ----
  dibujar <- function() {
    op <- par(mar = c(1, 1, 3, 1)); on.exit(par(op))
    M <- out$M
    leg <- if (M == 2) c("Régimen bajo", "Régimen alto") else paste("Régimen", seq_len(M))
    if (opaca) {
      # Caras opacas del color (claro) del régimen; persp dibuja de atrás hacia
      # adelante, así lo de adelante tapa lo de atrás y la red no se cruza.
      claro <- sapply(col_lin[-K], function(cc) {
        v <- grDevices::col2rgb(cc) / 255; v <- v + (1 - v) * 0.45
        grDevices::rgb(v[1], v[2], v[3])
      })
      caras <- matrix(rep(claro, each = N), nrow = K - 1, byrow = TRUE)
      persp(x = x, y = 0:N, z = t(Z), theta = theta, phi = phi, expand = 0.55,
            col = caras, border = adjustcolor("grey25", 0.6), lwd = 0.3,
            shade = 0.25, ltheta = -60, ticktype = "detailed", bg = fondo,
            xlab = "Trimestre del shock", ylab = "Horizonte (trimestres)", zlab = eje_z,
            main = titulo, cex.lab = 0.85, cex.axis = 0.7)
      legend("topleft", bty = "n", cex = 0.8, legend = leg, col = .colores_reg(M), pch = 15, pt.cex = 1.6)
      return(invisible(NULL))
    }
    # caras claras y bordes (la red) del color del régimen del trimestre
    caras <- matrix(rep(adjustcolor(col_lin[-K], 0.12), each = N), nrow = K - 1, byrow = TRUE)
    pm <- persp(x = x, y = 0:N, z = t(Z), theta = theta, phi = phi, expand = 0.55,
                col = caras, border = NA, ticktype = "detailed", bg = fondo,
                xlab = "Trimestre del shock", ylab = "Horizonte (trimestres)", zlab = eje_z,
                main = titulo, cex.lab = 0.85, cex.axis = 0.7, lwd = 0.4)
    # líneas a lo largo del horizonte (una por trimestre), coloreadas por régimen
    for (k in seq_len(K)) lines(trans3d(rep(x[k], N + 1), 0:N, Z[, k], pm), col = col_lin[k], lwd = 0.6)
    # líneas a lo largo del tiempo (una por horizonte)
    for (hh in 0:N) lines(trans3d(x, rep(hh, K), Z[hh + 1, ], pm), col = adjustcolor("grey15", 0.7), lwd = 0.6)
    legend("topleft", bty = "n", cex = 0.8, legend = leg, col = .colores_reg(M), lwd = 2)
  }
  dibujar()
  .abrir_pdf(paste0(archivo, ".pdf"), width = ancho, height = alto)
  dibujar()
  dev.off()

  # ---- HTML interactivo: superficie con malla ----
  tk <- .ticks_trim(x)
  vc <- matrix(rep(.valor_color(out)[idx], each = N + 1), nrow = N + 1)
  etq <- .etiqueta_trim(x)
  TX <- matrix(paste0(rep(etq, each = N + 1), "<br>Horizonte: ", rep(0:N, K),
                      "<br>", eje_z, ": ", sprintf("%.3f", as.vector(Z))), nrow = N + 1)
  lim_x <- range(x)
  g <- plot_ly(x = x, y = 0:N, z = Z, type = "surface", surfacecolor = vc, cmin = 0, cmax = 1,
               colorscale = .escala_plotly(out$M), colorbar = .barra_plotly(out$M),
               opacity = 0.85, text = TX, hoverinfo = "text",
               contours = list(
                 x = list(show = TRUE, start = floor(lim_x[1]), end = ceiling(lim_x[2]), size = 1,
                          color = "rgba(30,30,30,0.55)", width = 1),
                 y = list(show = TRUE, start = 0, end = N, size = 1,
                          color = "rgba(30,30,30,0.55)", width = 1))) %>%
    layout(title = titulo,
           scene = list(xaxis = list(title = "Trimestre del shock", tickvals = tk$at, ticktext = tk$lab),
                        yaxis = list(title = "Horizonte (trimestres)", dtick = 1),
                        zaxis = list(title = eje_z),
                        aspectratio = list(x = 2, y = 1, z = 0.8)))
  tryCatch(saveWidget(g, paste0(archivo, ".html"), selfcontained = FALSE),
           error = function(e) message("No se pudo guardar el HTML: ", conditionMessage(e)))
  g
}


# =====================================================================
# GIF: el ERPT a un horizonte fijo (por defecto h = 8) se va dibujando a lo
# largo del tiempo, trimestre a trimestre, con bandas de color según el régimen
# y una etiqueta con la fecha y el valor del trimestre actual.
# Requiere gifski. No usa 'res'.
# =====================================================================
guardar_gif_erpt_h <- function(out, h = 8, archivo = "gif_erpt_h8", cada = 1,
                               segundos_por_cuadro = 0.1, pausa_final = 3,
                               eje_y = "ERPT", titulo = NULL, intensidad = 0.30,
                               ancho = 1000, alto = 560, marcas = NULL, pausa_marca = 1) {
  if (!requireNamespace("gifski", quietly = TRUE)) {
    stop("Falta el paquete gifski. Instalalo una vez con: install.packages('gifski')")
  }
  out <- .normalizar_out(out)
  if (h > out$N) stop("El horizonte h no puede ser mayor que ", out$N, ".")
  M  <- out$M
  x  <- out$anio_dec
  z  <- out$erpt[h + 1, ]
  K  <- length(x)
  etq <- .etiqueta_trim(x)
  tk <- .ticks_trim(x)
  cols_banda <- sapply(seq_len(K), function(k) adjustcolor(.color_historia(out, k), intensidad))
  cols_punto <- sapply(seq_len(K), function(k) .color_historia(out, k))
  rango <- range(z, na.rm = TRUE); rango <- rango + c(-0.05, 0.12) * diff(rango)
  if (is.null(titulo)) titulo <- paste0(eje_y, " acumulado a h = ", h, " según el trimestre del shock")
  leg <- if (M == 2) c("Régimen bajo", "Régimen alto") else paste("Régimen", seq_len(M))

  # Marcas: trimestres elegidos ("3q1989", ...) que quedan etiquetados con fecha y valor
  im <- integer(0)
  if (!is.null(marcas)) {
    im <- match(marcas, etq)
    if (anyNA(im)) stop("No encontré: ", paste(marcas[is.na(im)], collapse = ", "),
                        ". Usá el formato '1q2002'.")
  }
  idx <- sort(unique(c(seq(1, K, by = cada), K, im)))
  n_pausa <- max(0, round(pausa_marca / segundos_por_cuadro))
  idx <- unlist(lapply(idx, function(k) if (k %in% im) rep(k, 1 + n_pausa) else k))   # pausa en cada marca
  idx <- c(idx, rep(K, max(0, round(pausa_final / segundos_por_cuadro))))           # pausa al final

  gifski::save_gif({
    for (k in idx) {
      par(mar = c(4.5, 4.5, 3, 1))
      plot(NA, xlim = range(x) + c(-0.125, 0.125), ylim = rango, xaxt = "n", xaxs = "i",
           xlab = "Trimestre del shock", ylab = paste0(eje_y, " (h = ", h, ")"), main = titulo)
      axis(1, at = tk$at, labels = tk$lab)
      rect(x[1:k] - 0.125, rango[1] - 1, x[1:k] + 0.125, rango[2] + 1,
           col = cols_banda[1:k], border = NA)
      abline(h = 0, col = "grey60")
      lines(x[1:k], z[1:k], col = "black", lwd = 1.8)
      points(x[k], z[k], pch = 19, cex = 1.6, col = cols_punto[k])
      # Etiqueta con fecha, régimen y valor del trimestre actual
      lab <- paste0(etq[k], "   ", .texto_estado(out, k), "   ", eje_y, " = ", sprintf("%.2f", z[k]))
      text(x[k], z[k], labels = sprintf("%.2f", z[k]), pos = 3, offset = 0.8, cex = 1.1, font = 2)
      mtext(lab, side = 3, line = 0.2, cex = 0.9, adj = 1)
      # Marcas ya alcanzadas: quedan fijas con su fecha y valor
      for (j in im[im <= k]) {
        segments(x[j], rango[1], x[j], z[j], col = "grey30", lty = 3)
        points(x[j], z[j], pch = 21, bg = cols_punto[j], col = "black", cex = 1.3)
        text(x[j], z[j], labels = paste0(etq[j], "\n", sprintf("%.2f", z[j])),
             pos = 3, offset = 0.9, cex = 0.85)
      }
      legend("topleft", bty = "n", cex = 0.85, legend = leg,
             col = adjustcolor(.colores_reg(M), intensidad + 0.2), pch = 15, pt.cex = 2)
    }
  }, gif_file = paste0(archivo, ".gif"), width = ancho, height = alto, delay = segundos_por_cuadro)
  message("GIF guardado en: ", normalizePath(paste0(archivo, ".gif")))
}


# =====================================================================
# GIF DE LA MANTA 3D (opaca, coloreada por régimen)
#   modo = "construir": la manta se va armando a lo largo del tiempo, trimestre a trimestre.
#   modo = "rotar":     la manta completa gira para verla desde distintos ángulos.
# 'marcas' (opcional): trimestres que quedan etiquetados con fecha y ERPT a h = h_marca.
# Requiere gifski. No usa 'res'.
# =====================================================================
guardar_gif_manta <- function(out, archivo = "gif_manta", modo = c("construir", "rotar"),
                              cada = 2, segundos_por_cuadro = 0.12, pausa_final = 3,
                              theta = -40, phi = 28, giro = 360, paso_giro = 5,
                              marcas = NULL, h_marca = out$N, eje_z = "ERPT", titulo = NULL,
                              ancho = 1000, alto = 750) {
  if (!requireNamespace("gifski", quietly = TRUE)) {
    stop("Falta el paquete gifski. Instalalo una vez con: install.packages('gifski')")
  }
  modo <- match.arg(modo)
  out <- .normalizar_out(out)
  N <- out$N; M <- out$M
  sel <- seq(1, ncol(out$erpt), by = cada)
  x   <- out$anio_dec[sel]
  Z   <- out$erpt[, sel, drop = FALSE]
  K   <- length(x)
  etq <- .etiqueta_trim(x)
  col_lin <- sapply(sel, function(k) .color_historia(out, k))
  claro <- sapply(col_lin[-K], function(cc) {
    v <- grDevices::col2rgb(cc) / 255; v <- v + (1 - v) * 0.45
    grDevices::rgb(v[1], v[2], v[3])
  })
  caras <- matrix(rep(claro, each = N), nrow = K - 1, byrow = TRUE)
  zlim <- range(Z, na.rm = TRUE)
  if (is.null(titulo)) titulo <- paste0(eje_z, " dependiente del estado")
  leg <- if (M == 2) c("Régimen bajo", "Régimen alto") else paste("Régimen", seq_len(M))

  im <- integer(0)
  if (!is.null(marcas)) {
    im <- match(marcas, .etiqueta_trim(out$anio_dec))
    if (anyNA(im)) stop("No encontré: ", paste(marcas[is.na(im)], collapse = ", "),
                        ". Usá el formato '1q2002'.")
    im <- sapply(im, function(k) which.min(abs(sel - k)))   # columna más cercana en la grilla
  }

  cuadro <- function(Zk, th, k_actual = NULL) {
    par(mar = c(1, 1, 3, 1))
    pm <- persp(x = x, y = 0:N, z = t(Zk), zlim = zlim, theta = th, phi = phi, expand = 0.55,
                col = caras, border = adjustcolor("grey25", 0.6), lwd = 0.3,
                shade = 0.25, ltheta = -60, ticktype = "detailed",
                xlab = "Trimestre del shock", ylab = "Horizonte (trimestres)", zlab = eje_z,
                main = titulo, cex.lab = 0.85, cex.axis = 0.7)
    legend("topleft", bty = "n", cex = 0.8, legend = leg, col = .colores_reg(M), pch = 15, pt.cex = 1.6)
    if (!is.null(k_actual)) mtext(etq[k_actual], side = 3, line = 0.2, adj = 1, cex = 1, font = 2)
    lim_k <- if (is.null(k_actual)) K else k_actual
    for (j in im[im <= lim_k]) {
      p <- trans3d(x[j], h_marca, Z[h_marca + 1, j], pm)
      points(p, pch = 21, bg = col_lin[j], col = "black", cex = 1.3)
      text(p, labels = paste0(etq[j], "\n", sprintf("%.2f", Z[h_marca + 1, j])),
           pos = 3, offset = 0.7, cex = 0.8, font = 2)
    }
  }

  n_fin <- max(0, round(pausa_final / segundos_por_cuadro))
  gifski::save_gif({
    if (modo == "construir") {
      for (k in c(2:K, rep(K, n_fin))) {
        Zk <- Z; if (k < K) Zk[, (k + 1):K] <- NA
        cuadro(Zk, theta, k)
      }
    } else {
      angulos <- theta + seq(0, giro, by = paso_giro)
      for (th in c(angulos, rep(tail(angulos, 1), n_fin))) cuadro(Z, th)
    }
  }, gif_file = paste0(archivo, ".gif"), width = ancho, height = alto, delay = segundos_por_cuadro)
  message("GIF guardado en: ", normalizePath(paste0(archivo, ".gif")))
}


# =====================================================================
# ACOTAR VALORES EXTREMOS SOLO PARA GRAFICAR
# Los gráficos (superficie 3D, manta, slider, GIF, red, 2D) muestran el ERPT de CADA
# trimestre, no un promedio ni una mediana: si en un trimestre la depreciación acumulada
# queda casi en cero, su ERPT se dispara y aplasta el resto del gráfico.
# Devuelve una copia de 'out' en la que cada valor se compara con los de los trimestres
# de SU MISMO RÉGIMEN (régimen de mayor peso) en ese horizonte: los que quedan fuera de
# [Q1 - k*IQR, Q3 + k*IQR] Y además se alejan de la mediana del régimen en más de
# 'umbral_abs' (en unidades de ERPT; 1 = 100 puntos de traslado) se llevan a ese límite.
# El umbral absoluto evita tocar diferencias chicas cuando los valores de un régimen están
# muy concentrados (p. ej. en los primeros horizontes): solo se acotan valores absurdos.
# NO cambia los resultados: las tablas (promedio y mediana) se calculan con el 'out' original.
# Uso: graficar_erpt_manta(acotar_para_graficar(out_...), archivo = "...")
# =====================================================================
acotar_para_graficar <- function(out, k = 5, umbral_abs = 1, min_trimestres = 8) {
  out0 <- out
  outn <- .normalizar_out(out)
  e <- out$erpt
  reg <- outn$regimen
  acotados <- list()
  for (h in seq_len(nrow(e))) {
    for (g in sort(unique(reg))) {
      idx <- which(reg == g)
      if (length(idx) < min_trimestres) next
      x  <- e[h, idx]
      q  <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
      ri <- q[2] - q[1]
      if (!is.finite(ri)) next
      med <- median(x, na.rm = TRUE)
      li <- min(q[1] - k * ri, med - umbral_abs)
      ls <- max(q[2] + k * ri, med + umbral_abs)
      fuera <- which(x < li | x > ls)
      if (length(fuera)) {
        acotados[[length(acotados) + 1]] <- data.frame(
          trimestre = colnames(e)[idx[fuera]], regimen = g, horizonte = h - 1,
          valor_original = round(x[fuera], 3),
          valor_graficado = round(pmin(pmax(x[fuera], li), ls), 3))
        e[h, idx[fuera]] <- pmin(pmax(x[fuera], li), ls)
      }
    }
  }
  if (length(acotados)) {
    tab <- do.call(rbind, acotados)
    tab <- tab[order(tab$trimestre, tab$horizonte), ]
    cat("Valores acotados SOLO para los gráficos (", nrow(tab), " de ", length(e),
        "; por régimen y horizonte: fuera de Q1 - ", k, "·IQR / Q3 + ", k,
        "·IQR y a más de ", umbral_abs, " de la mediana):\n", sep = "")
    print(tab, row.names = FALSE)
  } else {
    cat("No hay valores extremos para acotar (k = ", k, ", umbral_abs = ", umbral_abs, ").\n", sep = "")
  }
  out0$erpt <- e
  if (!is.null(out0$factor) && !is.null(out0$erpt_escalado)) out0$erpt_escalado <- e / out0$factor
  out0$acotado_para_graficar <- TRUE
  invisible(out0)
}
