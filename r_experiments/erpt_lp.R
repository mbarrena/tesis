# =====================================================================
# ERPT con proyecciones locales (lpirfs): relativo al shock, acumulado e incrementales,
# con bandas de 68% y 95% por bootstrap por bloques. Lineal o con umbral (2 regímenes).
#
# Uso (después de source("utils.r"), que define prepare_threshold_dummy):
#   source("erpt_lp.R")
#
#   # Lineal
#   lp_erpt <- erpt_lp(df_ERPT_Arg, endog = c("impp_usa", "E", "ipc", "pbird"),
#                      fixed_lags = 2, horizons = 8)
#   graficar_erpt_lp(lp_erpt, archivo = "erpt_lp_impp_usa_endog_p2")
#
#   # Con umbral: los mismos argumentos de umbral que run_lp_model
#   lp_erpt_u <- erpt_lp(df_ERPT_Arg, endog = c("impp_usa", "E", "ipc", "pbird"),
#                        fixed_lags = 2, horizons = 8,
#                        threshold_var = "ipc", upper_threshold = 0.80)
#   graficar_erpt_lp(lp_erpt_u, archivo = "erpt_lp_impp_usa_endog_p2_umbral80")
#
# Definiciones (en cada horizonte h; el shock es unitario: E sube 1 en el impacto):
#   relativo al shock = inflación acumulada hasta h / depreciación inicial
#   acumulado         = inflación acumulada hasta h / depreciación acumulada hasta h (transmisión)
#   persistencia      = depreciación acumulada / depreciación inicial
#   incremental       = ERPT en h menos ERPT en h-1 (en h = 0, el propio ERPT)
#
# Bandas: bootstrap por bloques móviles. Se remuestrean bloques de trimestres consecutivos
# (endógenas, exógenas y la variable de estado juntas), se reestiman las LP y se calculan
# las mismas medidas. Así las bandas valen también para los cocientes y para las respuestas
# acumuladas (sumar los límites de cada horizonte, como hace cumulative = TRUE en
# run_lp_model, no da un intervalo de confianza válido para la suma).
# =====================================================================

library(lpirfs)
library(ggplot2)

# Nombres de los regímenes de lp_nl con una dummy (use_logistic = FALSE):
# lpirfs llama s1 al régimen con dummy = 0 y s2 al régimen con dummy = 1.
.nombres_regimen <- c(s1 = "Régimen 1 (debajo del umbral)", s2 = "Régimen 2 (encima del umbral)")

# --- Una estimación de LP con las mismas opciones que run_lp_model (shock unitario) ---
.lp_estimar <- function(endog_data, exog_data, fixed_lags, horizons, lags_exog, trend,
                        newey_lags, switching = NULL, num_cores = NULL) {
  if (is.null(switching)) {
    lp_lin(
      endog_data,
      exog_data      = exog_data,
      lags_criterion = NaN,
      max_lags       = NA,
      lags_endog_lin = fixed_lags,
      trend          = trend,
      shock_type     = 1,
      confint        = 1.96,
      nw_lag         = newey_lags,
      nw_prewhite    = TRUE,
      hor            = horizons,
      lags_exog      = lags_exog,
      num_cores      = num_cores
    )
  } else {
    lp_nl(
      endog_data,
      exog_data      = exog_data,
      lags_criterion = NaN,
      max_lags       = NA,
      lags_endog_lin = fixed_lags,
      lags_endog_nl  = fixed_lags,
      trend          = trend,
      shock_type     = 1,
      confint        = 1.96,
      nw_lag         = newey_lags,
      nw_prewhite    = TRUE,
      hor            = horizons,
      lags_exog      = lags_exog,
      switching      = switching,
      use_logistic   = FALSE,
      use_hp         = FALSE,
      num_cores      = num_cores
    )
  }
}

# --- Las medidas de ERPT a partir de las respuestas de ipc y de E (horizontes 0..H) ---
.medidas_erpt <- function(resp_p, resp_e) {
  acum_p <- cumsum(resp_p)
  acum_e <- cumsum(resp_e)
  e0 <- resp_e[1]
  rel <- acum_p / e0
  acu <- acum_p / acum_e
  cbind(
    inflacion_acum          = acum_p,
    depreciacion_acum       = acum_e,
    erpt_relativo_shock     = rel,
    erpt_acumulado          = acu,
    persistencia            = acum_e / e0,
    erpt_relativo_shock_inc = c(rel[1], diff(rel)),
    erpt_acumulado_inc      = c(acu[1], diff(acu))
  )
}

# --- Medidas de cada régimen (o del modelo lineal) a partir de un resultado de lpirfs ---
.extraer <- function(res, i_p, i_e, umbral) {
  if (!umbral) {
    list(lineal = .medidas_erpt(res$irf_lin_mean[i_p, , i_e], res$irf_lin_mean[i_e, , i_e]))
  } else {
    list(s1 = .medidas_erpt(res$irf_s1_mean[i_p, , i_e], res$irf_s1_mean[i_e, , i_e]),
         s2 = .medidas_erpt(res$irf_s2_mean[i_p, , i_e], res$irf_s2_mean[i_e, , i_e]))
  }
}

# --- Índices de un bootstrap por bloques móviles ---
.indices_bloques <- function(n, bloque) {
  n_bloques <- ceiling(n / bloque)
  inicios <- sample.int(n - bloque + 1, n_bloques, replace = TRUE)
  idx <- unlist(lapply(inicios, function(s) s:(s + bloque - 1)))
  idx[seq_len(n)]
}

#' ERPT de una LP (lineal o con umbral) con bandas por bootstrap por bloques
#'
#' @param df data frame con las variables en diferencias logarítmicas (df_ERPT_Arg)
#' @param endog variables endógenas, en el orden de identificación (Cholesky)
#' @param exog variables exógenas (opcional)
#' @param fixed_lags cantidad de rezagos
#' @param horizons horizonte máximo (8 = dos años)
#' @param impulso, respuesta variable del shock ("E") y de precios ("ipc")
#' @param threshold_var, upper_threshold, lower_threshold, discard_threshold
#'        umbral, igual que en run_lp_model (NULL = modelo lineal)
#' @param reps réplicas del bootstrap (0 = sin bandas)
#' @param bloque largo de los bloques en trimestres (por defecto, horizons)
#' @param signifs niveles de las bandas
#' @param semilla semilla del bootstrap
erpt_lp <- function(df, endog, exog = NULL, fixed_lags, horizons = 8,
                    impulso = "E", respuesta = "ipc",
                    threshold_var = NULL, upper_threshold = 0.80, lower_threshold = NULL,
                    discard_threshold = NULL,
                    reps = 500, bloque = NULL, signifs = c(0.68, 0.95), semilla = 123,
                    lags_exog = NULL, trend = 0, newey_lags = NULL, mostrar = TRUE) {

  if (!(impulso %in% endog) || !(respuesta %in% endog)) {
    stop("'", impulso, "' y '", respuesta, "' tienen que estar entre las endógenas.")
  }
  i_e <- match(impulso, endog)
  i_p <- match(respuesta, endog)
  if (i_e > i_p) {
    warning("Con Cholesky, ", impulso, " va después de ", respuesta,
            ": el ipc no responde en el impacto. Revisá el orden.")
  }

  # ---- Umbral (misma lógica que run_lp_model) ----
  umbral <- !is.null(threshold_var)
  switching <- NULL
  if (umbral) {
    if (!exists("prepare_threshold_dummy")) stop("Falta source(\"utils.r\") antes de usar el umbral.")
    # Valor del umbral en unidades de la variable (el mismo cálculo que hace prepare_threshold_dummy)
    # (si se descartan datos por arriba del umbral superior, la dummy usa el inferior, y viceversa)
    p_umbral <- if (!is.null(discard_threshold) && discard_threshold == "upper") lower_threshold else
                if (!is.null(upper_threshold)) upper_threshold else lower_threshold
    valor_umbral <- unname(stats::quantile(df[[threshold_var]], p_umbral, na.rm = TRUE))
    pr <- prepare_threshold_dummy(df, threshold_var, upper_threshold, lower_threshold, discard_threshold)
    df <- pr$df
    switching <- as.numeric(pr$dummy)
    if (!all(switching %in% c(0, 1))) stop("La variable de estado tiene que ser 0/1.")
    cat("Umbral:", threshold_var, ">", round(valor_umbral, 4), "(percentil", round(100 * p_umbral, 1), ")\n")
    cat("Trimestres por régimen: debajo del umbral (dummy = 0) =", sum(switching == 0),
        "| encima del umbral (dummy = 1) =", sum(switching == 1), "\n")
  }

  endog_data <- data.frame(lapply(df[, endog, drop = FALSE], as.numeric))
  exog_data <- NULL
  if (!is.null(exog)) {
    exog_data <- data.frame(lapply(df[, exog, drop = FALSE], as.numeric))
    if (is.null(lags_exog)) lags_exog <- fixed_lags
  }
  completos <- stats::complete.cases(endog_data)
  if (!is.null(exog_data)) completos <- completos & stats::complete.cases(exog_data)
  if (!all(completos)) {
    stop("Hay valores faltantes en las variables del modelo (filas ",
         paste(head(which(!completos), 10), collapse = ", "), "...).")
  }

  # ---- Estimación puntual ----
  est <- .lp_estimar(endog_data, exog_data, fixed_lags, horizons, lags_exog, trend, newey_lags,
                     switching = switching)
  med <- .extraer(est, i_p, i_e, umbral)
  nombres <- if (umbral) .nombres_regimen[names(med)] else c(lineal = "Lineal")
  tabla <- do.call(rbind, lapply(names(med), function(r)
    data.frame(regimen = nombres[[r]], horizonte = 0:horizons, med[[r]], check.names = FALSE)))
  medidas <- colnames(med[[1]])

  # ---- Bootstrap por bloques ----
  if (reps > 0) {
    if (is.null(bloque)) bloque <- max(horizons, 4)
    n <- nrow(endog_data)
    set.seed(semilla)
    sims <- lapply(med, function(m) array(NA_real_, dim = c(reps, horizons + 1, length(medidas)),
                                          dimnames = list(NULL, NULL, medidas)))
    fallidas <- 0
    message("Bootstrap por bloques: ", reps, " réplicas, bloques de ", bloque, " trimestres...")
    pb <- utils::txtProgressBar(min = 0, max = reps, style = 3)
    for (rr in seq_len(reps)) {
      idx <- .indices_bloques(n, bloque)
      res_r <- tryCatch({
        e_r <- endog_data[idx, , drop = FALSE]; rownames(e_r) <- NULL
        x_r <- NULL
        if (!is.null(exog_data)) { x_r <- exog_data[idx, , drop = FALSE]; rownames(x_r) <- NULL }
        s_r <- if (umbral) switching[idx] else NULL
        b <- suppressMessages(suppressWarnings(
          .lp_estimar(e_r, x_r, fixed_lags, horizons, lags_exog, trend, newey_lags,
                      switching = s_r, num_cores = 1)))
        .extraer(b, i_p, i_e, umbral)
      }, error = function(e) NULL)
      if (is.null(res_r)) {
        fallidas <- fallidas + 1
      } else {
        for (r in names(sims)) sims[[r]][rr, , ] <- res_r[[r]]
      }
      utils::setTxtProgressBar(pb, rr)
    }
    close(pb)
    if (fallidas > 0) message(fallidas, " réplicas fallaron y se descartaron.")

    bandas <- do.call(rbind, lapply(names(sims), function(r) {
      d <- data.frame(regimen = nombres[[r]], horizonte = 0:horizons)
      for (s in signifs) {
        nivel <- round(100 * s)
        for (m in medidas) {
          d[[paste0(m, "_li", nivel)]] <- apply(sims[[r]][, , m], 2, stats::quantile, probs = (1 - s) / 2, na.rm = TRUE)
          d[[paste0(m, "_ls", nivel)]] <- apply(sims[[r]][, , m], 2, stats::quantile, probs = 1 - (1 - s) / 2, na.rm = TRUE)
        }
      }
      cerca <- colMeans(abs(sims[[r]][, , "depreciacion_acum"]) < 0.1, na.rm = TRUE)
      if (any(cerca > 0.05)) {
        message("Aviso (", nombres[[r]], "): en más del 5% de las réplicas la depreciación acumulada es casi cero en h = ",
                paste(which(cerca > 0.05) - 1, collapse = ", "),
                "; las bandas del ERPT acumulado en esos horizontes son poco confiables.")
      }
      d
    }))
    tabla <- cbind(tabla, bandas[, setdiff(names(bandas), c("regimen", "horizonte")), drop = FALSE])
  }

  attr(tabla, "impulso") <- impulso
  attr(tabla, "respuesta") <- respuesta
  attr(tabla, "signifs") <- signifs
  attr(tabla, "modelo") <- paste0("LP, ", fixed_lags, " rezagos",
                                   if (!is.null(exog)) paste0(", exógenas: ", paste(exog, collapse = ", ")) else "",
                                   if (umbral) paste0(", umbral: ", threshold_var, " > ",
                                                      formatC(valor_umbral, format = "f", digits = 3, decimal.mark = ","),
                                                      ", ", sum(switching == 1), " trimestres encima") else "")
  if (umbral) {
    attr(tabla, "umbral") <- c(percentil = p_umbral, valor = valor_umbral,
                               trimestres_debajo = sum(switching == 0), trimestres_encima = sum(switching == 1))
  }
  if (mostrar) {
    cat("ERPT de un shock unitario en", impulso, "sobre", respuesta, "(", attr(tabla, "modelo"), ")\n")
    for (r in unique(tabla$regimen)) {
      cat("\n", r, "\n", sep = "")
      print(round(tabla[tabla$regimen == r, c("horizonte", "erpt_relativo_shock", "erpt_acumulado", "persistencia",
                                               "erpt_relativo_shock_inc", "erpt_acumulado_inc")], 4), row.names = FALSE)
    }
  }
  invisible(tabla)
}

#' Gráfico de ERPT de las LP: relativo al shock y acumulado (y sus incrementales).
#' Con umbral, los dos regímenes van superpuestos en cada panel.
graficar_erpt_lp <- function(tabla, incremental = TRUE, titulo = NULL, archivo = NULL,
                             ancho = 10, alto = NULL, niveles = NULL) {
  # niveles: qué bandas dibujar, entre las calculadas en erpt_lp (p. ej. niveles = 0.68).
  #          NULL = todas. No hace falta volver a estimar.
  paneles <- c(erpt_relativo_shock = "Relativo al shock",
               erpt_acumulado = "Acumulado (transmisión)")
  if (incremental) {
    paneles <- c(paneles,
                 erpt_relativo_shock_inc = "Relativo al shock: incremental",
                 erpt_acumulado_inc = "Acumulado (transmisión): incremental")
  }
  signifs <- attr(tabla, "signifs")
  if (!is.null(niveles)) {
    faltan <- setdiff(round(niveles, 2), round(signifs, 2))
    if (length(faltan) > 0) stop("Las bandas de ", paste(faltan, collapse = ", "),
                                 " no están calculadas; las disponibles son: ", paste(signifs, collapse = ", "),
                                 ". Para otras, volvé a correr erpt_lp con signifs = c(...).")
    signifs <- signifs[round(signifs, 2) %in% round(niveles, 2)]
  }
  signifs <- sort(signifs, decreasing = TRUE)   # la más ancha primero (abajo), la más angosta encima
  # Intensidad de cada banda: la más angosta, más oscura
  etiquetas <- paste0(round(100 * signifs), "%")
  intensidad <- if (length(signifs) == 1) "oscura" else c(rep("clara", length(signifs) - 1), "oscura")
  texto_bandas <- paste0("Bandas: ", paste0(rev(etiquetas),
                         if (length(signifs) > 1) paste0(" (más ", rev(intensidad), ")") else "",
                         collapse = " y "), ", por bootstrap por bloques")
  regimenes <- unique(tabla$regimen)
  varios <- length(regimenes) > 1

  largo <- do.call(rbind, lapply(names(paneles), function(m) {
    d <- data.frame(panel = paneles[[m]], regimen = tabla$regimen,
                    horizonte = tabla$horizonte, valor = tabla[[m]])
    for (s in signifs) {
      nivel <- round(100 * s)
      li <- paste0(m, "_li", nivel)
      if (li %in% names(tabla)) {
        d[[paste0("li", nivel)]] <- tabla[[li]]
        d[[paste0("ls", nivel)]] <- tabla[[paste0(m, "_ls", nivel)]]
      }
    }
    d
  }))
  largo$panel <- factor(largo$panel, levels = unname(paneles))
  largo$regimen <- factor(largo$regimen, levels = regimenes)

  if (is.null(titulo)) {
    titulo <- paste0("ERPT: shock en ", attr(tabla, "impulso"), " sobre ", attr(tabla, "respuesta"),
                     " (", attr(tabla, "modelo"), ")")
  }

  if (!varios) {
    # Un solo modelo: bandas azules (68% oscura, 95% clara), línea negra
    colores <- setNames(ifelse(intensidad == "oscura", "#7aa6cc", "#c6d8ea"), etiquetas)
    g <- ggplot(largo, aes(x = horizonte, y = valor))
    for (s in signifs) {
      nivel <- as.character(round(100 * s))
      li <- paste0("li", nivel); ls <- paste0("ls", nivel)
      if (li %in% names(largo)) {
        g <- g + geom_ribbon(aes(ymin = .data[[li]], ymax = .data[[ls]], fill = !!paste0(nivel, "%")))
      }
    }
    g <- g +
      geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4) +
      geom_line(aes(colour = "Estimación"), linewidth = 1) +
      scale_fill_manual(values = colores, breaks = rev(etiquetas), name = NULL) +
      scale_colour_manual(values = c("Estimación" = "black"), name = NULL) +
      labs(caption = texto_bandas)
  } else {
    # Dos regímenes superpuestos: un color por régimen; banda de 68% más opaca que la de 95%
    col_reg <- setNames(c("steelblue", "firebrick")[seq_along(regimenes)], regimenes)
    alfas <- setNames(ifelse(intensidad == "oscura", 0.28, 0.12), as.character(round(100 * signifs)))
    g <- ggplot(largo, aes(x = horizonte, y = valor, colour = regimen, fill = regimen))
    for (s in signifs) {
      nivel <- as.character(round(100 * s))
      li <- paste0("li", nivel); ls <- paste0("ls", nivel)
      if (li %in% names(largo)) {
        g <- g + geom_ribbon(aes(ymin = .data[[li]], ymax = .data[[ls]]), colour = NA, alpha = alfas[[nivel]])
      }
    }
    g <- g +
      geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4) +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = col_reg, name = NULL) +
      scale_fill_manual(values = col_reg, name = NULL) +
      labs(caption = texto_bandas)
  }

  g <- g +
    facet_wrap(~ panel, ncol = 2, scales = "free_y") +
    scale_x_continuous(breaks = sort(unique(tabla$horizonte))) +
    labs(title = titulo, x = "Horizonte (trimestres)", y = "ERPT") +
    theme_bw(base_size = 12) +
    theme(legend.position = "top", plot.title = element_text(hjust = 0.5, size = 13),
          strip.background = element_blank(), strip.text = element_text(size = 12))

  if (is.null(alto)) alto <- if (incremental) 8 else 4.5
  options(repr.plot.width = ancho, repr.plot.height = alto)
  print(g)
  if (!is.null(archivo)) ggsave(paste0(archivo, ".pdf"), g, width = ancho, height = alto)
  invisible(NULL)
}
