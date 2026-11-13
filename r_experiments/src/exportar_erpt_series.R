# =====================================================================
# exportar_erpt_series.R
# Exporta a Excel el ERPT de cada período (trimestre o año) a los horizontes
# elegidos: relativo a la depreciación inicial, transmisión y persistencia,
# con un gráfico de líneas por modelo (PNG insertado en la hoja y guardado aparte).
#
# Requiere erpt_historia_v38.R cargado antes (usa .normalizar_out).
# Uso: ver las celdas de ejemplo al final del archivo.
# =====================================================================

if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")

# Busca un objeto out_... en memoria; si no está, lo lee de <carpeta>/<nombre>.rds
cargar_out <- function(nombre, carpeta = ".") {
  if (exists(nombre, envir = .GlobalEnv)) return(get(nombre, envir = .GlobalEnv))
  f <- file.path(carpeta, paste0(nombre, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  stop("No encuentro '", nombre, "' ni en memoria ni como ", f,
       ". Corré la celda de erpt_historia de ese modelo (con load_presaved = TRUE) ",
       "o indicá la carpeta donde está el .rds.", call. = FALSE)
}

# Arma la tabla período por período para un modelo
.serie_erpt <- function(rel, acum, horizontes) {
  rel  <- .normalizar_out(rel)
  acum <- .normalizar_out(acum)
  if (!identical(colnames(rel$erpt), colnames(acum$erpt)))
    stop("El objeto relativo y el acumulado no tienen los mismos períodos: ¿son del mismo modelo?",
         call. = FALSE)
  if (identical(rel$definicion, "acumulado") || identical(acum$definicion, "relativo"))
    stop("Parece que 'relativo' y 'transmision' están invertidos.", call. = FALSE)
  N <- nrow(rel$erpt) - 1
  if (any(horizontes > N)) stop("Pediste un horizonte mayor que N = ", N, ".", call. = FALSE)

  M     <- ncol(rel$pesos)
  fecha <- if (!is.null(rel$anio_dec)) rel$anio_dec else seq_len(ncol(rel$erpt))
  df <- data.frame(periodo = colnames(rel$erpt), fecha = round(fecha, 2),
                   regimen = rel$regimen, stringsAsFactors = FALSE)
  df[[if (M == 2) "peso_reg_alto" else paste0("peso_reg_", M)]] <- round(rel$pesos[, M], 4)
  for (h in horizontes) df[[paste0("relativo_h", h)]]    <- unname(rel$erpt[h + 1, ])
  for (h in horizontes) df[[paste0("transmision_h", h)]] <- unname(acum$erpt[h + 1, ])
  for (h in horizontes) df[[paste0("persistencia_h", h)]] <-
    df[[paste0("relativo_h", h)]] / df[[paste0("transmision_h", h)]]
  attr(df, "M")     <- M
  attr(df, "anual") <- identical(rel$frecuencia, "anual")
  df
}

# Gráfico de líneas: arriba relativo, abajo transmisión; sombreado = régimen alto
.graficar_serie <- function(df, horizontes, titulo, archivo_png, unidad) {
  M     <- attr(df, "M")
  medio <- if (attr(df, "anual")) 0.5 else 0.125
  cols  <- c("#1f4e9c", "#d95f02", "#1b9e77", "#7570b3")[seq_along(horizontes)]
  png(archivo_png, width = 11, height = 7.5, units = "in", res = 150)
  op <- par(mfrow = c(2, 1), mar = c(3, 4.5, 2.5, 1), oma = c(0, 0, 3, 0))
  for (medida in c("relativo", "transmision")) {
    Y <- do.call(cbind, lapply(horizontes, function(h) df[[paste0(medida, "_h", h)]]))
    # Eje acotado a los percentiles 2-98 (solo el gráfico; la tabla tiene todos los valores)
    lim <- quantile(Y, c(0.02, 0.98), na.rm = TRUE)
    lim <- range(lim + c(-0.1, 0.1) * diff(lim), 0, 1)
    plot(NA, xlim = range(df$fecha), ylim = lim, xlab = "", ylab = "ERPT", las = 1,
         main = if (medida == "relativo") "Relativo a la depreciación inicial"
                else "Transmisión (inflación acumulada / depreciación acumulada)")
    for (k in seq_len(M)[-1]) {
      sel <- df$regimen == k
      if (any(sel)) rect(df$fecha[sel] - medio, lim[1] - 10, df$fecha[sel] + medio, lim[2] + 10,
                         col = adjustcolor("firebrick", alpha.f = 0.10 + 0.08 * (k - 2)), border = NA)
    }
    abline(h = 0, col = "grey40"); abline(h = 1, lty = 3, col = "grey40")
    for (i in seq_along(horizontes)) {
      y     <- Y[, i]
      fuera <- which(y < lim[1] | y > lim[2])
      lines(df$fecha, pmin(pmax(y, lim[1]), lim[2]), col = cols[i], lwd = 1.6)
      if (length(fuera)) points(df$fecha[fuera], ifelse(y[fuera] > lim[2], lim[2], lim[1]),
                                pch = 4, col = cols[i])
    }
    legend("topleft", legend = paste0("h = ", horizontes, " ", unidad), col = cols, lwd = 2,
           bty = "n", horiz = TRUE, cex = 0.85)
  }
  mtext(titulo, outer = TRUE, line = 1.4, cex = 1.1, font = 2)
  mtext("Sombreado: régimen de alta inflación.  ×: valor fuera de escala (ver la tabla).",
        outer = TRUE, line = 0.2, cex = 0.8)
  par(op)
  dev.off()
}

# modelos: lista con nombre; cada elemento = list(relativo = out, transmision = out_acum)
exportar_erpt_series <- function(modelos, horizontes, archivo, carpeta_graficos = NULL) {
  wb    <- openxlsx::createWorkbook()
  est_h <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#D9E1F2",
                                 wrapText = TRUE, border = "bottom")
  est_n <- openxlsx::createStyle(numFmt = "0.000")

  openxlsx::addWorksheet(wb, "Leeme")
  openxlsx::writeData(wb, "Leeme", data.frame(Notas = c(
    "Una hoja por modelo; una fila por período (trimestre o año en que ocurre el shock).",
    "relativo_hX = inflación acumulada hasta X / depreciación inicial.",
    "transmision_hX = inflación acumulada hasta X / depreciación acumulada hasta X.",
    "persistencia_hX = relativo / transmisión = depreciación acumulada / depreciación inicial (< 0,8: la depreciación se revierte y la transmisión pierde sentido).",
    "regimen = régimen de mayor peso en ese período; peso_reg_alto = peso del régimen alto (alfa).",
    "Valores en unidades originales, los mismos que grafican las mantas, los GIF y los sliders (erpt_historia).",
    "Gráfico: eje Y acotado a los percentiles 2-98 solo para el dibujo; la tabla tiene todos los valores.",
    paste("Generado:", format(Sys.time(), "%Y-%m-%d %H:%M")))))
  openxlsx::setColWidths(wb, "Leeme", cols = 1, widths = 130)

  dir_png <- if (is.null(carpeta_graficos)) tempdir() else carpeta_graficos
  if (!dir.exists(dir_png)) dir.create(dir_png, recursive = TRUE)

  salida <- list()
  for (nm in names(modelos)) {
    m  <- modelos[[nm]]
    df <- .serie_erpt(m$relativo, m$transmision, horizontes)
    unidad <- if (attr(df, "anual")) "años" else "trimestres"
    hoja   <- substr(gsub("[^A-Za-z0-9_ .=()-]", "", nm), 1, 31)

    openxlsx::addWorksheet(wb, hoja)
    openxlsx::writeData(wb, hoja, df, headerStyle = est_h)
    cols_num <- grep("^(relativo|transmision|persistencia)_h", names(df))
    openxlsx::addStyle(wb, hoja, est_n, rows = 2:(nrow(df) + 1), cols = cols_num, gridExpand = TRUE)
    openxlsx::setColWidths(wb, hoja, cols = seq_along(df), widths = 13)
    openxlsx::freezePane(wb, hoja, firstRow = TRUE)

    png_f <- file.path(dir_png, paste0("serie_erpt_", gsub("[^A-Za-z0-9]+", "_", nm), ".png"))
    .graficar_serie(df, horizontes, nm, png_f, unidad)
    openxlsx::insertImage(wb, hoja, png_f, width = 11, height = 7.5, units = "in",
                          startRow = 2, startCol = ncol(df) + 2)
    cat(sprintf("%-32s %3d períodos | gráfico: %s\n", nm, nrow(df), png_f))
    salida[[nm]] <- df
  }
  openxlsx::saveWorkbook(wb, archivo, overwrite = TRUE)
  cat("\nListo:", archivo, "\n")
  invisible(salida)
}

# =====================================================================
# CELDAS DE EJEMPLO (copiar a la notebook; ajustar nombres si difieren)
# =====================================================================
#
# --- Trimestral (notebook de 2 regímenes, 1950-2023) ---
# source("src/erpt_historia_v38.R")
# source("src/exportar_erpt_series.R")
# modelos_trim <- list(
#   "TVAR p2 (principal)" = list(relativo    = cargar_out("out_impp_usa_tvar_p2"),
#                                transmision = cargar_out("out_impp_usa_tvar_p2_acum")),
#   "TVAR p3 (robustez)"  = list(relativo    = cargar_out("out_impp_usa_tvar_p3"),
#                                transmision = cargar_out("out_impp_usa_tvar_p3_acum")))
# series_trim <- exportar_erpt_series(modelos_trim, horizontes = c(4, 8),
#                                     archivo = "ERPT_series_trimestrales.xlsx",
#                                     carpeta_graficos = "GRAFICOS_TVAR_VLSTAR_Trimestrales")
#
# --- Anual (notebook de Ferreres; lee de los .rds las salidas de fuentes combinadas) ---
# source("src/exportar_erpt_series.R")
# modelos_anual <- list(
#   "Comb TVAR p1 (principal)"  = list(relativo    = cargar_out("out_impp_usa_tvar_p1_anual_comb"),
#                                      transmision = cargar_out("out_impp_usa_tvar_p1_acum_anual_comb")),
#   "Ferr TVAR p1 (principal)"  = list(relativo    = cargar_out("out_impp_usa_tvar_p1_anual_ferr"),
#                                      transmision = cargar_out("out_impp_usa_tvar_p1_acum_anual_ferr")),
#   "Ferr TVAR-t p1 (robustez)" = list(relativo    = cargar_out("out_impp_usa_tvar_t_p1"),
#                                      transmision = cargar_out("out_impp_usa_tvar_t_p1_acum")))
# series_anual <- exportar_erpt_series(modelos_anual, horizontes = c(1, 2, 4),
#                                      archivo = "ERPT_series_anuales.xlsx",
#                                      carpeta_graficos = "GRAFICOS_TVAR_VLSTAR_Anuales")
