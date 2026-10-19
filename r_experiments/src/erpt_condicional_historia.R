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
                          regimen_alto = 2,     # régimen de alta inflación (el de mayor pass-through)
                          definicion = c("impacto", "acumulado"),
                          modelo     = NULL,    # si no lo encuentra en 'res', pasalo a mano
                          escalas    = NULL,    # desvíos usados para escalar (se buscan en res$...$scaling_sd)
                          semilla    = 123,
                          ncores     = max(1, parallel::detectCores() - 1)) {

  definicion <- match.arg(definicion)

  # ---- Modelo estructural y escalas ----
  if (is.null(modelo)) {
    modelo <- .buscar(res, function(x) inherits(x, "stvar") &&
                        !identical(x$model$identification, "reduced_form"))
  }
  if (is.null(modelo)) stop("No encontré el modelo estructural en 'res'. Pasalo con modelo = ...")
  if (is.null(escalas)) escalas <- .buscar_nombre(res, "scaling_sd")
  if (is.null(escalas)) escalas <- sapply(data[, endog], sd, na.rm = TRUE)  # misma cuenta que fit_generic_sstvar

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
  etiqueta <- paste0(data$año, "T", q)

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
  alpha <- modelo$transition_weights[t_shock - p, regimen_alto]

  # ---- GIRF condicional a cada historia ----
  # Mismo shock (1 desvío estructural) y misma semilla en todas las historias:
  # la variación entre trimestres refleja solo el estado de la economía.
  message("Estimando ", K, " GIRFs (una por trimestre), ", R1, " simulaciones cada una...")
  g <- GIRF(modelo,
            which_shocks     = shock_num,
            shock_size       = 1,
            N                = N,
            R1               = R1,
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

  # ---- Chequeo: historias dominadas por cada régimen ----
  dom_alto <- alpha > 0.9
  dom_bajo <- alpha < 0.1
  chequeo <- data.frame(
    horizonte          = 0:N,
    reg_bajo_escalado  = rowMeans(erpt_escalado[, dom_bajo, drop = FALSE]),
    reg_alto_escalado  = rowMeans(erpt_escalado[, dom_alto, drop = FALSE]),
    reg_bajo_original  = rowMeans(erpt[, dom_bajo, drop = FALSE]),
    reg_alto_original  = rowMeans(erpt[, dom_alto, drop = FALSE])
  )
  cat("\nTrimestres con alpha > 0.9:", sum(dom_alto), "| con alpha < 0.1:", sum(dom_bajo), "\n")
  cat("Factor de conversión a unidades originales (sd_", resp_var, "/ sd_", shock_var, "): ",
      round(sd_p / sd_e, 4), "\n\n", sep = "")
  print(round(chequeo, 4))

  invisible(list(erpt = erpt, erpt_escalado = erpt_escalado, alpha = alpha,
                 anio_dec = anio_dec[t_shock], etiqueta = etiqueta[t_shock],
                 N = N, definicion = definicion, chequeo = chequeo,
                 factor = sd_p / sd_e))
}

# =====================================================================
# Gráficos
# =====================================================================
graficar_erpt <- function(out, archivo = "erpt", horizontes_2d = c(4, out$N)) {
  N <- out$N
  write.csv(out$erpt, paste0(archivo, "_matriz.csv"))
  write.csv(data.frame(trimestre = out$etiqueta, alpha_reg_alto = out$alpha),
            paste0(archivo, "_alpha.csv"), row.names = FALSE)

  # ---- 3D interactivo (color = peso del régimen de alta inflación) ----
  color_alpha <- matrix(rep(out$alpha, each = N + 1), nrow = N + 1)
  g3d <- plot_ly(x = out$anio_dec, y = 0:N, z = out$erpt, type = "surface",
                 surfacecolor = color_alpha, cmin = 0, cmax = 1,
                 colorscale = list(c(0, "steelblue"), c(1, "firebrick")),
                 colorbar = list(title = "\u03b1 régimen alto")) %>%
    layout(title = "ERPT dependiente del estado",
           scene = list(xaxis = list(title = "Trimestre del shock"),
                        yaxis = list(title = "Horizonte (trimestres)"),
                        zaxis = list(title = "ERPT")))
  saveWidget(g3d, paste0(archivo, "_3d.html"), selfcontained = TRUE)

  # ---- 3D estático para la tesis ----
  pdf(paste0(archivo, "_3d.pdf"), width = 9, height = 7)
  persp(x = out$anio_dec, y = 0:N, z = t(out$erpt),
        theta = -40, phi = 25, expand = 0.6, ticktype = "detailed",
        col = "lightsteelblue", border = NA, shade = 0.4,
        xlab = "Trimestre del shock", ylab = "Horizonte (trimestres)", zlab = "ERPT",
        main = "ERPT dependiente del estado")
  dev.off()

  # ---- 2D complementario: ERPT a horizontes fijos + alpha sombreado ----
  pdf(paste0(archivo, "_2d.pdf"), width = 10, height = 5)
  par(mar = c(4, 4, 3, 4))
  rango <- range(out$erpt[horizontes_2d + 1, ], na.rm = TRUE)
  plot(NA, xlim = range(out$anio_dec), ylim = rango,
       xlab = "Trimestre del shock", ylab = "ERPT", main = "ERPT a horizontes seleccionados")
  polygon(c(out$anio_dec, rev(out$anio_dec)),
          c(rango[1] + out$alpha * diff(rango), rep(rango[1], length(out$alpha))),
          col = adjustcolor("firebrick", 0.15), border = NA)
  cols <- c("steelblue", "black", "darkgreen", "purple")
  for (i in seq_along(horizontes_2d)) {
    lines(out$anio_dec, out$erpt[horizontes_2d[i] + 1, ], col = cols[i], lwd = 1.5)
  }
  axis(4, at = rango[1] + c(0, 0.5, 1) * diff(rango), labels = c(0, 0.5, 1))
  mtext("\u03b1 régimen alto (sombreado)", side = 4, line = 2.5)
  legend("topleft", bty = "n", lwd = 1.5, col = cols[seq_along(horizontes_2d)],
         legend = paste0("h = ", horizontes_2d, " trimestres"))
  dev.off()

  g3d
}
