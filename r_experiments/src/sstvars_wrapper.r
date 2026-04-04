library(tsDyn)
library(sstvars)

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
  estim_method = "two-phase"
) {
  Y <- data[, endog]
  seeds <- 1:nrounds
  thr_params <- c(which(colnames(Y) == thr_var), thr_lag)

  tvar_model <- fitSTVAR(
    data = Y,
    p = feature_lag,
    M = n_regimes,
    weight_function = weight_function, # threshold
    weightfun_pars = thr_params,
    nrounds = nrounds,
    seeds = seeds,
    estim_method = estim_method
  )

  cat("========= Unstructured fitted model ==========\n")
  summary(tvar_model)

  fit_struct <- fitSSTVAR(
    stvar = tvar_model,
    identification = "heteroskedasticity"
  )

  tryCatch(
    {
      cat("Saving trained models to models/ folder\n")
      saveRDS(tvar_model, paste0("models/", filename))
      saveRDS(fit_struct, paste0("models/", "structured", "_", filename))
    },
    error = function(e) {
      cat("!!! No se pudo guardar modelos pre-estimados. GUARDAR MANUALMENTE corriendo save_trained_models(models, models/", filename, "). Error: ", e$message, "\n")
    }
  )
  cat("========= Structured fitted model (heteroskedasticity) ==========\n")
  tryCatch(
    {
      summary(fit_struct)
    },
    error = function(e) {
      cat("!!! No se pudo obtener summary del modelo estructurado. Error: ", e$message, "\n")
    }
  )

  return(c(tvar_model, fit_struct))
}

load_trained_model <- function(filename, structured = TRUE) {
  return(readRDS(paste0("models/", ifelse(structured, paste0("structured", "_"), ""), filename)))
}

save_trained_models <- function(models, filename) {
  saveRDS(models[[1]], paste0("models/", filename))
  saveRDS(models[[2]], paste0("models/", "structured", "_", filename))
  cat("Modelos guardados exitosamente en carpeta models/.")
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
  load_presaved = FALSE
) {
  # validate_input_cols(data, endog, shock_var, response_vars)

  tvar_model <- NULL
  fit_struct <- NULL
  if (load_presaved) {
    tryCatch(
      {
        tvar_model <- load_trained_model(filename, structured = FALSE)
        fit_struct <- load_trained_model(filename, structured = TRUE)
        cat("Se cargó exitosamente modelos preguardados ", filename, " de la carpeta models\n")
      },
      error = function(e) {
        cat("!!! No se pudo cargar modelo pre-estimado. Correr de nuevo con load_presaved = FALSE. Error: ", e$message, "\n")
        return(NULL)
      }
    )
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
      nrounds = nrounds
    )
    tvar_model <- fit_models[1]
    fit_struct <- fit_models[2]
  }

  shock_var_idx <- which(colnames(fit_struct$data) == shock_var)
  response_vars_idx <- which(colnames(fit_struct$data) %in% response_vars)

  scale_matrix <- sapply(response_vars_idx, function(x) {
    c(shock_var_idx, x, irf_shocks)
  })

  if (cumulative) {
    cumulative_idx <- which(response_vars %in% response_vars)
  } else {
    cumulative_idx <- NULL
  }

  if (is.null(irf_horizons)) {
    if (periodicity == "trimestral") {
      irf_horizons <- 8
    } else if (periodicity == "anual") {
      irf_horizons <- 2
    }
  }

  if (run_tests) {
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
    cat("\n\n")
  }

  if (diagnostic_plots) {
    cat(" >>> DIAGNOSTIC PLOTS\n")
    plot(fit_struct)
    diagnostic_plot(tvar_model, type = "series", resid_type = "standardized", maxlag = irf_horizons)
  }

  irf_response <- list(regimes = list())

  for (i in 1:n_regimes) {
    cat("+++++++++ RÉGIMEN ", i, "+++++++++\n")
    cat(" >>> GIRF\n")
    girf_erpt <- GIRF(
      fit_struct,
      N = irf_horizons, # horizonte
      which_shocks = shock_var_idx, # shock al Tipo de Cambio (variable 2)
      which_cumulative = cumulative_idx, # acumula la respuesta de inflación (variable 3)
      init_regime = i,
      scale = scale_matrix,
      scale_type = "instant",
      R1 = nrounds
    )
    plot(girf_erpt)

    cat(" >>> Linear IRF (1 shock)")
    irf_struct <- linear_IRF(
      fit_struct,
      N = irf_horizons, # horizonte de 8 periodos
      regime = i,
      which_cumulative = cumulative_idx, # acumula las respuestas de la inflación (variable 3)
      scale = scale_matrix,
      ncores = 2
    )
    plot(irf_struct)

    reg_irfs <- list(girf = girf_erpt, lirf = irf_struct)
    irf_response$regimes[[i]] <- reg_irfs
  }

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
  load_presaved = FALSE
) {
  weight_function <- "threshold"
  filename <- paste0(
    "threshold_tvar_",
    thr_var, "_",
    thr_lag, "_",
    nrounds,
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
    load_presaved = load_presaved
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
  load_presaved = FALSE
) {
  weight_function <- "logistic"
  filename <- paste0(
    "vlstar_",
    thr_var, "_",
    thr_lag, "_",
    nrounds,
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
    load_presaved = load_presaved
  )
}

print_irf <- function(res) {
  print(res$regimes[[1]]$lirf)
}

print_girf <- function(res) {
  print(res$regimes[[1]]$girf)
}
