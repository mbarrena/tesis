cat(">>> sstvars_wrapper v5 cargado (ind_Student con identificación por no normalidad)\n")
library(tsDyn)
library(sstvars)
library(digest)

#' Run Local Projection Model
#'
#' This function estimates a local projection (LP) model using the `lpirfs` package,
#' allowing for endogenous and exogenous variables, different lag structures,
#' and cumulative impulse responses.
#'
#' @param data A data frame containing the time series data.
#' @param endog A character vector specifying the names of endogenous variables.
#' @param feature_lag An integer indicating a fixed number of lags for the model features.
#' @param threshold_var A character vector specifying the name of the threshold variable's column
#' @param threshold_var An integer indicating a fixed number of lags for the transition variable
#' @param filename A character string indicating the name of the file to save the trained model to. Saves to /models folder
#' @param nrounds An integer indicating the number of rounds for the estimation. Default 500
#' @param seeds (optional) An integer indicating the number of rounds for the estimation. Default 1:nrounds
#' @param estim_method (optional) A character indicating the estimation method. Default "two-phase"
#'
#' @return A list containing impulse response function (IRF) results.
#'         The function also generates plots and prints results.
#'
#' @details
#' - The function runs the `lp_lin()` function from the `lpirfs` package to estimate a local projection model.
#' - If `cumulative = TRUE`, the cumulative impulse responses are computed.
#' - The function generates plots using `plotlp()` and formats results using `pretty_results()`.
#'
#' @examples
#' \dontrun{
#' data <- your_dataframe
#' results <- fit_generic_stvar()
#' }
#'
#' @export
fit_generic_sstvar <- function(
  data,
  endog,
  feature_lag,
  n_regimes,
  weight_function,
  thr_var,
  thr_lag,
  filename,
  nrounds,
  seeds = NULL,
  estim_method = "two-phase",
  cond_dist = "Gaussian",
  min_obs_coef = 3
) {
  Y <- data[, endog]

  # --- Escalado por desvío estándar (sobre los datos ya cappeados) ---
  Y_sd <- sapply(Y, sd, na.rm = TRUE)
  Y <- sweep(Y, 2, Y_sd, FUN = "/")
  cat("Desvíos usados para escalar:\n")
  print(Y_sd)
  # --------------------------------------------------------------------

  seeds <- if (is.null(seeds)) 1:nrounds else seeds
  estim_method <- if (is.null(estim_method)) "two-phase" else estim_method
  thr_idx <- which(colnames(Y) == thr_var)
  # mlogit pide una lista (variables de transición y número de rezagos);
  # threshold y logistic piden c(variable, rezago)
  thr_params <- if (identical(weight_function, "mlogit")) {
    list(vars = thr_idx, lags = thr_lag)
  } else {
    c(thr_idx, thr_lag)
  }

  args_fit <- list(
    data = Y,
    p = feature_lag,
    M = n_regimes,
    weight_function = weight_function, # threshold
    weightfun_pars = thr_params,
    cond_dist = cond_dist,
    nrounds = nrounds,
    seeds = seeds,
    estim_method = estim_method,
    ncores = 6
  )
  # min_obs_coef solo actúa en la estimación en tres etapas; se pasa solo en ese caso
  # (así la estimación en dos etapas queda exactamente igual que en las versiones anteriores)
  if (identical(estim_method, "three-phase")) args_fit$min_obs_coef <- min_obs_coef
  tvar_model <- do.call(fitSTVAR, args_fit)


  cat("========= Unstructured fitted model ==========\n")
  summary(tvar_model)

  if (cond_dist %in% c("ind_Student", "ind_skewed_t")) {
    # Con shocks t independientes el modelo ya queda identificado por no normalidad
    # (Virolainen): no se usa la identificación por heteroscedasticidad.
    fit_struct <- tryCatch(
      fitSSTVAR(stvar = tvar_model, identification = "non-Gaussianity"),
      error = function(e) {
        message("fitSSTVAR(non-Gaussianity) no corrió (", conditionMessage(e),
                "); se usa el modelo estimado, que ya está identificado por no normalidad.")
        m <- tvar_model
        m$model$identification <- "non-Gaussianity"
        m
      })
    cat("========= Structured fitted model (non-Gaussianity) ==========\n")
  } else {
    fit_struct <- fitSSTVAR(
      stvar = tvar_model,
      identification = "heteroskedasticity"
    )
    cat("========= Structured fitted model (heteroskedasticity) ==========\n")
  }

  summary(fit_struct)

  tryCatch(
    {
      cat("Saving trained models to models/ folder\n")
      saveRDS(tvar_model, paste0("models/", filename))
      saveRDS(fit_struct, paste0("models/", "structured", "_", filename))
    },
    error = function(e) {
      cat("!!! No se pudo guardar modelos pre-estimados. GUARDAR MANUALMENTE corriendo save_trained_models(models, models/", filename, "). Error: ", e$message, "\n")
      return(NULL)
    }
  )
  return(list(tvar_model = tvar_model, fit_struct = fit_struct))
}

load_trained_model <- function(filename, structured = TRUE) {
  return(readRDS(paste0("models/", ifelse(structured, paste0("structured", "_"), ""), filename)))
}

save_trained_models <- function(models, filename) {
  saveRDS(models[[1]], paste0("models/", filename))
  saveRDS(models[[2]], paste0("models/", "structured", "_", filename))
  cat("Modelos guardados exitosamente en carpeta models/.")
}

run_tests_sstvar <- function(models, response_vars, cumulative, irf_horizons = NULL, irf_shocks = 1, nrounds = 500) {
  tvar_model <- models$tvar_model
  fit_struct <- models$fit_struct

  if (cumulative) {
    cumulative_idx <- which(colnames(fit_struct$data) %in% response_vars)
  } else {
    cumulative_idx <- NULL
  }

  if (is.null(irf_horizons)) {
    if (periodicity == "trimestral") {
      irf_horizons <- 8
    } else if (periodicity == "anual") {
      irf_horizons <- 2
    } else if (periodicity == "mensual") {
      irf_horizons <- 24
    }
  }

  cat(" >>> Portmanteau\n")
  port_res <- Portmanteau_test(tvar_model) # autocorrelación
  print(port_res)
  cat("\n\n")
  # LR_test(fit_struct) # lineal vs no lineal
  cat(" >>> Historical Decomposition\n")
  decomp_res <- hist_decomp(fit_struct)
  print(decomp_res)
  cat("\n\n")
  cat(" >>> GFEVD\n")
  gefvd_res <- suppressMessages(GFEVD(
    fit_struct,
    N = irf_horizons,
    shock_size = irf_shocks,
    which_cumulative = cumulative_idx,
    R1 = nrounds
  ))
  print(gefvd_res)
}

run_sstvar_fit_tests_irf <- function(
  data,
  endog,
  weight_function = weight_function,
  filename = filename,
  periodicity, # c('trimestral','anual')
  feature_lag,
  n_regimes,
  cumulative,
  thr_var, # string
  thr_lag,
  shock_var, # string
  response_vars, # list(strings)
  run_tests = TRUE, # run tests
  irf_horizons = NULL, # Max periods after to look to in IRF (default 2 years)
  irf_shocks = 1, # Shocks unitarios irf (default 1)
  diagnostic_plots = FALSE, # Print diagnostic plots?
  nrounds = 500,
  load_presaved = FALSE,
  cond_dist = "Gaussian",
  estim_method = "two-phase",
  min_obs_coef = 3
) {
  # validate_input_cols(data, endog, shock_var, response_vars)

  # Desvíos usados para escalar Y = data[, endog] antes de estimar (mismos que calcula
  # fit_generic_sstvar internamente). Se recalculan acá para que estén disponibles
  # incluso cuando load_presaved = TRUE, y para no depender del return de fit_generic_sstvar.
  scaling_sd <- sapply(data[, endog], sd, na.rm = TRUE)

  tvar_model <- NULL
  fit_struct <- NULL
  if (load_presaved) {
    tvar_model <- load_trained_model(filename, structured = FALSE)
    fit_struct <- load_trained_model(filename, structured = TRUE)
    cat("Se cargó exitosamente modelos preguardados ", filename, " de la carpeta models\n")
  } else {
    fit_models <- fit_generic_sstvar(
      data,
      endog,
      feature_lag,
      n_regimes,
      weight_function,
      thr_var,
      thr_lag,
      filename,
      nrounds = nrounds,
      estim_method = estim_method,
      cond_dist = cond_dist,
      min_obs_coef = min_obs_coef
    )
    tvar_model <- fit_models$tvar_model
    fit_struct <- fit_models$fit_struct
  }

  if (run_tests) {
    cat("\n >>> TESTS\n")
    cat("PARA CORRER LOS TESTS CORRER LA SIGUIENTE FUNCION:\n")
    cat("run_tests_sstvar(res, c('", response_vars, "'),", cumulative, ",", ifelse(is.null(irf_horizons), "NULL", irf_horizons), ",", irf_shocks, ",", nrounds, ")\n")
  }

  shock_var_idx <- which(colnames(fit_struct$data) == shock_var)
  response_vars_idx <- which(colnames(fit_struct$data) %in% response_vars)

  scale_matrix <- sapply(response_vars_idx, function(x) {
    c(shock_var_idx, x, irf_shocks)
  })

  if (cumulative) {
    cumulative_idx <- response_vars_idx
  } else {
    cumulative_idx <- NULL
  }

  if (is.null(irf_horizons)) {
    if (periodicity == "trimestral") {
      irf_horizons <- 8
    } else if (periodicity == "anual") {
      irf_horizons <- 2
    } else if (periodicity == "mensual") {
      irf_horizons <- 24
    }
  }

  if (diagnostic_plots) {
    cat("\n >>> DIAGNOSTIC PLOTS\n")
    plot(fit_struct)
    diagnostic_plot(tvar_model, type = "series", resid_type = "standardized", maxlag = irf_horizons)
    plot.new()
  }

  irf_response <- list(regimes = list())

  for (i in 1:n_regimes) {
    cat("+++++++++ RÉGIMEN ", i, "+++++++++\n")
    cat(" >>> GIRF (", irf_shocks, " shock units)\n")
    girf_erpt <- GIRF(
      fit_struct,
      N = irf_horizons, # horizonte
      which_shocks = shock_var_idx, # shock al Tipo de Cambio (variable 2)
      which_cumulative = cumulative_idx, # acumula la respuesta de inflación (variable 3)
      init_regime = i,
      scale = scale_matrix,
      scale_type = "instant",
      R1 = nrounds,
      ncores = 6
    )
    plot(girf_erpt)
    plot.new()

    cat(" >>> Linear IRF (", irf_shocks, " shock units)\n")
    irf_struct <- linear_IRF(
      fit_struct,
      N = irf_horizons, # horizonte de 8 periodos
      regime = i,
      which_cumulative = cumulative_idx, # acumula las respuestas de la inflación (variable 3)
      scale = scale_matrix,
      ncores = 6
    )
    plot(irf_struct)
    plot.new()

    reg_irfs <- list(girf = girf_erpt, lirf = irf_struct)
    irf_response$regimes[[i]] <- reg_irfs
  }

  irf_response$scaling_sd <- scaling_sd

  # Guardar también los modelos estimados (necesario para explorar óptimos alternativos)
  irf_response$tvar_model <- tvar_model   # modelo sin identificar (contiene todas las rondas)
  irf_response$fit_struct <- fit_struct   # modelo identificado por heteroscedasticidad
  irf_response$filename   <- filename

  return(irf_response)
}

run_threshold_tvar <- function(
  data,
  endog,
  periodicity, # c('trimestral','anual')
  feature_lag,
  n_regimes,
  cumulative,
  thr_var, # string
  thr_lag,
  shock_var, # string
  response_vars, # list(strings)
  run_tests = TRUE, # run tests
  irf_horizons = NULL, # Max periods after to look to in IRF (default 2 years)
  irf_shocks = 1, # Shocks unitarios irf (default 1)
  diagnostic_plots = FALSE, # Print diagnostic plots?
  nrounds = 500,
  load_presaved = FALSE,
  cond_dist = "Gaussian",   # "Gaussian" (default) o "Student" (errores t de Student)
  estim_method = "two-phase",  # "two-phase" (default) o "three-phase"
  min_obs_coef = 3             # tres etapas: mínimo de observaciones por régimen (ver run_vlstar)
) {
  weight_function <- "threshold"
  include_params <- c("data","endog","periodicity","feature_lag","n_regimes","thr_var","thr_lag","nrounds")
  # cond_dist entra al hash solo si no es gaussiano: así los modelos gaussianos ya guardados
  # conservan su nombre de archivo y se siguen cargando con load_presaved = TRUE
  if (cond_dist != "Gaussian") include_params <- c(include_params, "cond_dist")
  # Igual que cond_dist: solo entran al hash si no son los valores por defecto,
  # para que los modelos ya guardados conserven su nombre de archivo
  if (estim_method != "two-phase") include_params <- c(include_params, "estim_method", "min_obs_coef")
  params_to_hash <- as.list(environment())
  params_to_hash <- params_to_hash[names(params_to_hash) %in% include_params]

  filename <- paste0(
    "threshold_tvar_",
    if (cond_dist != "Gaussian") paste0(tolower(cond_dist), "_") else "",
    if (estim_method != "two-phase") paste0("3fases_moc", min_obs_coef, "_") else "",
    feature_lag, "_",
    n_regimes, "_",
    nrounds, "_",
    digest(params_to_hash),
    ".rds"
  )

  run_sstvar_fit_tests_irf(
    data = data,
    endog = endog,
    weight_function = weight_function,
    filename = filename,
    periodicity = periodicity,
    feature_lag = feature_lag,
    n_regimes = n_regimes,
    cumulative = cumulative,
    thr_var = thr_var,
    thr_lag = thr_lag,
    shock_var = shock_var,
    response_vars = response_vars,
    run_tests = run_tests,
    irf_horizons = irf_horizons,
    irf_shocks = irf_shocks,
    diagnostic_plots = diagnostic_plots,
    nrounds = nrounds,
    load_presaved = load_presaved,
    cond_dist = cond_dist,
    estim_method = estim_method,
    min_obs_coef = min_obs_coef
  )
}


run_vlstar <- function(
  data,
  endog,
  periodicity, # c('trimestral','anual')
  feature_lag,
  n_regimes,
  cumulative,
  thr_var, # string
  thr_lag,
  shock_var, # string
  response_vars, # list(strings)
  run_tests = TRUE, # run tests
  irf_horizons = NULL, # Max periods after to look to in IRF (default 2 years)
  irf_shocks = 1, # Shocks unitarios irf (default 1)
  diagnostic_plots = FALSE, # Print diagnostic plots?
  nrounds = 500,
  load_presaved = FALSE,
  cond_dist = "Gaussian",   # "Gaussian" (default) o "Student" (errores t de Student)
  estim_method = "two-phase",  # "two-phase" (default) o "three-phase"
  min_obs_coef = 3             # tres etapas: mínimo de observaciones por régimen (ver run_vlstar)
) {
  weight_function <- "logistic"
  include_params <- c("data","endog","periodicity","feature_lag","n_regimes","thr_var","thr_lag","nrounds")
  # cond_dist entra al hash solo si no es gaussiano: así los modelos gaussianos ya guardados
  # conservan su nombre de archivo y se siguen cargando con load_presaved = TRUE
  if (cond_dist != "Gaussian") include_params <- c(include_params, "cond_dist")
  # Igual que cond_dist: solo entran al hash si no son los valores por defecto,
  # para que los modelos ya guardados conserven su nombre de archivo
  if (estim_method != "two-phase") include_params <- c(include_params, "estim_method", "min_obs_coef")
  params_to_hash <- as.list(environment())
  params_to_hash <- params_to_hash[names(params_to_hash) %in% include_params]

  filename <- paste0(
    "vlstar_",
    if (cond_dist != "Gaussian") paste0(tolower(cond_dist), "_") else "",
    if (estim_method != "two-phase") paste0("3fases_moc", min_obs_coef, "_") else "",
    feature_lag, "_",
    n_regimes, "_",
    nrounds, "_",
    digest(params_to_hash),
    ".rds"
  )

  run_sstvar_fit_tests_irf(
    data = data,
    endog = endog,
    weight_function = weight_function,
    filename = filename,
    periodicity = periodicity,
    feature_lag = feature_lag,
    n_regimes = n_regimes,
    cumulative = cumulative,
    thr_var = thr_var,
    thr_lag = thr_lag,
    shock_var = shock_var,
    response_vars = response_vars,
    run_tests = run_tests,
    irf_horizons = irf_horizons,
    irf_shocks = irf_shocks,
    diagnostic_plots = diagnostic_plots,
    nrounds = nrounds,
    load_presaved = load_presaved,
    cond_dist = cond_dist,
    estim_method = estim_method,
    min_obs_coef = min_obs_coef
  )
}

print_irf <- function(res) {
  for (i in 1:length(res$regimes)) {
    cat("+++++++++ RÉGIMEN ", i, "+++++++++\n")
    print(res$regimes[[i]]$lirf)
  }
}

print_girf <- function(res) {
  for (i in 1:length(res$regimes)) {
    cat("+++++++++ RÉGIMEN ", i, "+++++++++\n")
    print(res$regimes[[i]]$girf)
  }
}


#' Imprime el pass-through normalizado (respuesta / tamaño del shock inicial),
#' para todos los horizontes disponibles, tanto de GIRF como de IRF lineal.
#'
#' @param res El objeto devuelto por run_threshold_tvar() o run_vlstar()
#' @param endog El mismo vector de variables endógenas que le pasaste a esa función,
#'              en el mismo orden (ej. c("Psoja_USA","E","ipc","pbird"))
#' @param shock_var Nombre de la variable que recibe el shock (default "E")
#' @param response_var Nombre de la variable de respuesta a normalizar (default "ipc")
print_passthrough <- function(res, endog, shock_var = "E", response_var = "ipc") {
  shock_idx    <- which(endog == shock_var)
  response_idx <- which(endog == response_var)

  if (length(shock_idx) == 0 || length(response_idx) == 0) {
    stop("shock_var o response_var no están en el vector endog que pasaste")
  }

  if (is.null(res$scaling_sd)) {
    stop("res$scaling_sd no está presente. Volvé a correr run_threshold_tvar()/run_vlstar() con la versión actualizada de sstvars_wrapper.r")
  }

  # --- Desescalar: volver a las unidades originales antes de normalizar ---
  sd_shock    <- res$scaling_sd[shock_var]
  sd_response <- res$scaling_sd[response_var]
  factor_desescalado <- sd_response / sd_shock
  # --------------------------------------------------------------------

  for (i in seq_along(res$regimes)) {
    cat("+++++++++ RÉGIMEN ", i, "+++++++++\n")

    # --- GIRF ---
    girf_obj <- res$regimes[[i]]$girf
    girf_pe  <- girf_obj$girf_res[[1]]$point_est   # matriz: horizonte x variable
    shock_size_girf <- girf_pe[1, shock_idx]        # fila 1 = horizonte 0 (impacto)
    girf_response    <- girf_pe[, response_idx]
    girf_normalized  <- (girf_response / shock_size_girf) * factor_desescalado

    cat("\n>>> GIRF normalizada (", response_var, "/ shock inicial de", shock_var, ") - unidades originales\n")
    print(data.frame(
      horizonte    = 0:(length(girf_normalized) - 1),
      pass_through = round(girf_normalized, 4)
    ))

    # --- Linear IRF ---
    lirf_obj <- res$regimes[[i]]$lirf
    lirf_pe  <- lirf_obj$point_est                  # array [variable, shock, horizonte]
    shock_size_lirf <- lirf_pe[shock_idx, shock_idx, 1]
    lirf_response    <- lirf_pe[response_idx, shock_idx, ]
    lirf_normalized  <- (lirf_response / shock_size_lirf) * factor_desescalado

    cat("\n>>> IRF lineal normalizada (", response_var, "/ shock inicial de", shock_var, ") - unidades originales\n")
    print(data.frame(
      horizonte    = 0:(length(lirf_normalized) - 1),
      pass_through = round(lirf_normalized, 4)
    ))
    cat("\n")
  }
}


# =====================================================================
# STVAR con función de transición LOGIT MULTINOMIAL ("mlogit")
# Transición suave que admite 3 o más regímenes (la logística de run_vlstar solo admite 2).
# El peso de cada régimen depende de thr_var rezagada thr_lag períodos:
#   alpha_m,t ∝ exp(gamma_m' z_{t-1}),  z_{t-1} = (1, thr_var_{t-1}),  gamma_M = 0.
# No tiene un umbral c ni un gamma únicos: los puntos donde cambia el régimen dominante
# se calculan con resumen_modelo() de erpt_historia_v27.R.
# =====================================================================
run_mlogit <- function(
  data,
  endog,
  periodicity,
  feature_lag,
  n_regimes,
  cumulative,
  thr_var,
  thr_lag,
  shock_var,
  response_vars,
  run_tests = TRUE,
  irf_horizons = NULL,
  irf_shocks = 1,
  diagnostic_plots = FALSE,
  nrounds = 500,
  load_presaved = FALSE,
  cond_dist = "Gaussian",
  estim_method = "two-phase",
  min_obs_coef = 3
) {
  weight_function <- "mlogit"
  include_params <- c("data","endog","periodicity","feature_lag","n_regimes","thr_var","thr_lag","nrounds")
  if (cond_dist != "Gaussian") include_params <- c(include_params, "cond_dist")
  # Igual que cond_dist: solo entran al hash si no son los valores por defecto,
  # para que los modelos ya guardados conserven su nombre de archivo
  if (estim_method != "two-phase") include_params <- c(include_params, "estim_method", "min_obs_coef")
  params_to_hash <- as.list(environment())
  params_to_hash <- params_to_hash[names(params_to_hash) %in% include_params]

  filename <- paste0(
    "mlogit_",
    if (cond_dist != "Gaussian") paste0(tolower(cond_dist), "_") else "",
    if (estim_method != "two-phase") paste0("3fases_moc", min_obs_coef, "_") else "",
    feature_lag, "_",
    n_regimes, "_",
    nrounds, "_",
    digest(params_to_hash),
    ".rds"
  )

  run_sstvar_fit_tests_irf(
    data = data,
    endog = endog,
    weight_function = weight_function,
    filename = filename,
    periodicity = periodicity,
    feature_lag = feature_lag,
    n_regimes = n_regimes,
    cumulative = cumulative,
    thr_var = thr_var,
    thr_lag = thr_lag,
    shock_var = shock_var,
    response_vars = response_vars,
    run_tests = run_tests,
    irf_horizons = irf_horizons,
    irf_shocks = irf_shocks,
    diagnostic_plots = diagnostic_plots,
    nrounds = nrounds,
    load_presaved = load_presaved,
    cond_dist = cond_dist,
    estim_method = estim_method,
    min_obs_coef = min_obs_coef
  )
}
