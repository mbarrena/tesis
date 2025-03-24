library(lpirfs)
library(gridExtra)
library(ggplot2)
library(readxl)
library(httr)
library(IRdisplay)

library(reshape2)
library(purrr) 

#' Run Local Projection Model
#'
#' This function estimates a local projection (LP) model using the `lpirfs` package, 
#' allowing for endogenous and exogenous variables, different lag structures, 
#' and cumulative impulse responses.
#'
#' @param data A data frame containing the time series data.
#' @param endog A character vector specifying the names of endogenous variables.
#' @param max_lags An integer indicating the maximum number of lags for lags_crterion.
#' @param signif A numeric value (0.95 or 0.68) specifying the confidence level.
#' @param lags_criterion (optional) A string indicating the lag length criterion ('AICc', 'AIC' or 'BIC') (default: 'AIC')
#' @param exog (optional) A character vector specifying the names of exogenous variables (default: NULL).
#' @param newey_lags (optional) An integer specifying the Newey-West lag selection (default: NULL).
#' @param horizons (optional) An integer defining the forecast horizons (default: 10).
#' @param lags_exog (optional) An integer indicating the number of lags for exogenous variables (default: NULL, which sets it to match `max_lags`).
#' @param trend (optional) An integer indicating the trend types to include: 0 = no trend, 1 = include linear trend, 2 = include linear and quadratic trend (default: 0).
#' @param cumulative (optional) A logical value indicating whether to compute cumulative impulse responses (default: FALSE).
#' @param threshold_var (optional) A character vector specifying the name of the threshold column. If it has only 0 and 1 values, the threshold variable can be decomposed via the Hodrick-Prescott filter. Otherwise, it is decomposed via a logistic function. (default: NULL)
#'
#' @return A list containing impulse response function (IRF) results from the `lp_lin()` function.
#'         The function also generates plots and prints formatted results.
#'
#' @details 
#' - The function runs the `lp_lin()` function from the `lpirfs` package to estimate a local projection model.
#' - If `cumulative = TRUE`, the cumulative impulse responses are computed.
#' - The function generates plots using `plotlp()` and formats results using `pretty_results()`.
#'
#' @examples
#' \dontrun{
#' data <- your_dataframe
#' results <- run_lp_model(
#'   data = data,
#'   endog = c("GDP", "Inflation"),
#'   exog = c("InterestRate"),
#'   max_lags = 3,
#'   lags_criterion = 'AIC',
#'   newey_lags = 4,
#'   horizons = 10,
#'   signif = 0.95,
#'   lags_exog = 2,
#'   trend = 2,
#'   cumulative = FALSE
#' )
#' }
#'
#' @export
run_lp_model <- function(data, endog, exog=NULL, max_lags, lags_criterion = 'AIC', newey_lags = NULL, horizons = 10, signif, lags_exog = NULL, trend = 0, cumulative = FALSE, threshold_var = NULL) {
  # Map confidence levels to the corresponding values
  confint_map <- c("0.95" = 1.96, "0.68" = 1)

  # Check if the provided significance level is valid
  if (!as.character(signif) %in% names(confint_map)) {
    stop("Invalid significance level. Use 0.05 for 95% or 0.32 for 68%.")
  }

  if (!is.null(exog)) {
    exog_data <- data[, exog, drop = FALSE]
    if (is.null(lags_exog)) {
        lags_exog = max_lags
    }
  } else {
    exog_data = NULL
  }

  # Select endogenous variables
  endog_data <- data[, endog, drop = FALSE]
  
  # Convert to numeric (ensure proper model fitting)
  endog_data <- data.frame(lapply(endog_data, as.numeric))
  results_lin <- NULL
  unique_count <- NULL

  if(is_null(threshold_var)) {
    # Run the local projections model
    results_lin <- lp_lin(
      endog_data, 
      exog_data      = exog_data,
      lags_criterion = 'AIC',
      max_lags = max_lags,
      lags_endog_lin = NA,
      trend          = trend,  
      shock_type     = 1,  
      confint        = confint_map[as.character(signif)], 
      nw_lag         = newey_lags,
      nw_prewhite    = TRUE, 
      hor            = horizons,
      lags_exog = lags_exog
    )
  } else {
    switching <- data[, threshold_var, drop = TRUE]
    use_logistic <- TRUE
    if (all(switching %in% c(0, 1))) {
      display_html("Threshold var is binary: will not use logistic decomposition")
      use_logistic <- FALSE
    } else {
      display_html("Threshold var is not binary: using logistic decomposition")
    }
    unique_count <- length(unique(switching))

    results_lin <- lp_nl(
      endog_data, 
      exog_data      = exog_data,
      lags_criterion = 'AIC',
      max_lags = max_lags,
      lags_endog_lin = NA,
      lags_endog_nl = NA,
      trend          = trend,  
      shock_type     = 1,  
      confint        = confint_map[as.character(signif)], 
      nw_lag         = newey_lags,
      nw_prewhite    = TRUE, 
      hor            = horizons,
      lags_exog = lags_exog,
      switching = switching,
      use_logistic = use_logistic,
      use_hp = FALSE
    )
  }

  title_text <- paste0(
    "LocalProjection ", 
    if (!is.null(exog)) "(with exog)" else "(without exog)", 
    if (cumulative) " - Cumulative" else "",
    " - signif ", signif
  )

  if (cumulative) {
    for (i in 1:dim(results_lin$irf_lin_mean)[3]) {
      # Cumsum mean IRFs
      results_lin$irf_lin_mean[,,i] <- t(apply(results_lin$irf_lin_mean[,,i], 1, cumsum))
      
      # Cumsum error bands to match cumulative mean
      results_lin$irf_lin_low[,,i] <- t(apply(results_lin$irf_lin_low[,,i], 1, cumsum))
      results_lin$irf_lin_up[,,i]  <- t(apply(results_lin$irf_lin_up[,,i], 1, cumsum))
    }
  }

  df_specs_summary <- specs_summary(results_lin)
  display_html("Run configurations:")
  display(df_specs_summary)

  df_chosen_lags_list <- lag_summary(results_lin, endog)

  for (shock in names(df_chosen_lags_list)) {
    display_html(paste0("Lags for ", shock))  # Print title
    display(df_chosen_lags_list[[shock]])  # Print corresponding dataframe
  }

  if (is_null(threshold_var)) {
    plotlp_lin(results_lin, endog = endog, title_text = title_text)
    # Process and print results with cumulative option
    pretty_results_lin(results_lin = results_lin, endog_vars = endog)
  } else {
    plotlp_nl(results_lin, endog = endog, title_text = title_text, unique_count = unique_count, horizon_count = horizons+1)
    # Process and print results with cumulative option
    #pretty_results_nl(results_lin = results_lin, endog_vars = endog, unique_count)
  }

  return(results_lin)
}

plotlp_lin <- function(results_lin, endog, title_text) {
  # Generate plots
  linear_plots <- plot_lin(results_lin)

  display_html(title_text)
  
  # Show all plots
  lin_plots_all <- sapply(linear_plots, ggplotGrob)

  final_plot <- marrangeGrob(lin_plots_all, nrow = length(endog), ncol = length(endog), 
                             top = grid::textGrob(title_text, gp = grid::gpar(fontsize = 14, fontface = "bold")))

  print(final_plot)
  return(final_plot)
}

plotlp_nl <- function(results_lin, endog, title_text, unique_count, horizon_count) {
  # Generate plots
  linear_plots <- plot_nl(results_lin)
  display_html(title_text)
  
  # Flatten the linear_plots object manually using lapply to avoid unlist() confusion
  lin_plots_all <- lapply(1:unique_count, function(j) {
    # Extract all horizons for each gg_sN
    plot_list <- linear_plots[[paste0("gg_s", j)]]
    
    # Ensure that each gg_sN has horizons, e.g., gg_s1h1, gg_s1h2, etc.
    print(paste("Structure of gg_s", j, ":", str(plot_list)))
    
    # Return the list of plots for each group (e.g., gg_s1, gg_s2, ...)
    return(plot_list)
  })
  
  # Create a list to hold the plots for each impulse-response pair
  final_plots <- list()
  
  # Arrange the plots for each impulse-response pair side by side
  for (i in 1:length(endog)) {
    # Collect the ith plot from each gg_sN for all horizons
    plots_for_pair <- lapply(1:unique_count, function(j) {
      # Loop over horizons (h1, h2, ...) for each gg_sN
      plots_for_horizon <- lin_plots_all[[j]]
      
      # Ensure that the index i corresponds to the correct horizon (e.g., i = 1, 2, ...)
      plot <- plots_for_horizon[[i]]  # Extract the ith plot for this group
      
      return(plot)  # Return the plot for the ith index
    })
    
    # Flatten the list of plots for each impulse-response pair to avoid nested lists
    plots_for_pair <- unlist(plots_for_pair, recursive = TRUE)
    
    # Dynamically set ncol to the unique_count for side-by-side plots
    final_plots[[i]] <- do.call(marrangeGrob, c(list(grobs = plots_for_pair), 
                                               list(nrow = 1, ncol = unique_count, 
                                                    top = grid::textGrob(paste(title_text, "Impulse-Response Pair", i), 
                                                                         gp = grid::gpar(fontsize = 14, fontface = "bold")))))
  }
  
  # Print each final plot (side-by-side layout for each impulse-response pair)
  for (plot in final_plots) {
    print(plot)
  }
  
  return(final_plots)
}


pretty_results_lin <- function(results_lin, endog_vars) {
  # Extract IRF values, lower and upper bounds
  irf_array <- results_lin$irf_lin_mean
  lower_bound_array <- results_lin$irf_lin_low
  upper_bound_array <- results_lin$irf_lin_up

  # Get dimensions
  dims <- dim(irf_array)
  n_endog <- dims[1]  # Number of endogenous variables (responses)
  n_horizons <- dims[2]  # Number of horizons
  n_impulses <- dims[3]  # Number of impulses

  # Ensure the array is flattened correctly
  irf_values <- as.vector(aperm(irf_array, c(3, 1, 2)))  # Reorder to match impulse-response-horizon
  lower_values <- as.vector(aperm(lower_bound_array, c(3, 1, 2)))
  upper_values <- as.vector(aperm(upper_bound_array, c(3, 1, 2)))

  # Create index grid with correctly ordered impulse, response, and horizon
  indices <- expand.grid(
    impulse = 1:n_impulses,  # third dimension
    response = 1:n_endog,    # first dimension
    horizon = 1:n_horizons   # second dimension
  )

  # Create the data frame with corrected order
  irf_df <- data.frame(
    indices,
    irf_value = irf_values,
    lower_bound = lower_values,
    upper_bound = upper_values
  )


  # Exclude rows where impulse == response
  irf_df <- subset(irf_df, impulse != response)

  # Order by impulse, response, then horizon
  irf_df <- irf_df[order(irf_df$impulse, irf_df$response, irf_df$horizon), ]

  # Remove row names (drop index numbers)
  row.names(irf_df) <- NULL

  irf_df <- as.data.frame(irf_df)

  # Map impulse and response numbers to names
  irf_named_df <- irf_df
  irf_named_df$impulse <- endog_vars[irf_named_df$impulse]
  irf_named_df$response <- endog_vars[irf_named_df$response]

  # Split into a list of dataframes by impulse-response pair
  irf_list <- split(irf_named_df, list(irf_named_df$impulse, irf_named_df$response), drop = TRUE)

  # Print all dataframes in the list
  lapply(irf_list, display)
  return(irf_named_df)
}

lag_summary <- function(res, endog_vars) {
  # Get all shocks (names of the elements in chosen_lags)
  shock_names <- names(res$specs$chosen_lags)

  # Iterate over each shock and create a dataframe
  df_list <- lapply(shock_names, function(shock) {
    chosen_lags_list <- res$specs$chosen_lags[[shock]]  # Extract the list of matrices

    # Convert matrices to a dataframe
    df_chosen_lags <- do.call(cbind, chosen_lags_list)
    df_chosen_lags <- as.data.frame(df_chosen_lags)  # Ensure it's a dataframe

    # Rename columns to indicate horizons
    colnames(df_chosen_lags) <- paste0("H_", seq_along(chosen_lags_list))

    # Add endogenous variable names
    df_chosen_lags$Variable <- endog_vars  

    # Reorder to have 'Variable' as the first column
    df_chosen_lags <- df_chosen_lags[, c("Variable", colnames(df_chosen_lags)[1:(ncol(df_chosen_lags) - 1)])]

    return(df_chosen_lags)
  })

  # Assign shock names to the list elements
  names(df_list) <- shock_names

  return(df_list)
}

specs_summary <- function(res) {
  specs <- res$specs  # Extract specs
  
  # Create a one-row dataframe
  df_specs <- data.frame(
    lags_endog_lin = if (!is.null(specs$lags_endog_lin)) specs$lags_endog_lin else NA,
    lags_criterion = specs$lags_criterion,
    max_lags = specs$max_lags,
    trend = specs$trend,
    shock_type = specs$shock_type,
    confint = names(specs$confint),  # Extract confidence level name
    hor = specs$hor,
    use_nw = specs$use_nw,
    nw_prewhite = specs$nw_prewhite,
    nw_lag = ifelse(is.null(specs$nw_lag), NA, specs$nw_lag),
    adjust_se = specs$adjust_se,
    use_twosls = specs$use_twosls,
    model_type = specs$model_type,
    starts = specs$starts,
    ends = specs$ends,
    column_names = paste(specs$column_names, collapse = ", "),  # Convert vector to a single string
    endog = specs$endog,
    stringsAsFactors = FALSE  # Prevent automatic factor conversion
  )
  
  return(df_specs)
}

makeLogColumns <- function(lista, data) {
  df <- data  # Make a copy of the data frame
  for (c in lista) {
    df[[c]] <- log(as.numeric(df[[c]]))  # Convert to numeric and take log
  }
  return(df)
}

makeDiffColumns <- function(lista, data, factor = 1, drop = TRUE) {
  df <- data  # Make a copy of the data frame
  for (c in lista) {
    df[[c]] <- c(rep(NA, factor), diff(df[[c]], differences = factor))  # Compute differences
  }
  if (drop) {
    df <- df[(factor + 1):nrow(df), , drop = FALSE]  # Drop first `factor` rows
  }
  return(df)
}

renameColumnOfDataframe <- function(dataframe, oldColumnName, newColumnName) {
  df <- dataframe  # Make a copy of the dataframe
  colnames(df)[colnames(df) == oldColumnName] <- newColumnName  # Rename column
  return(df)
}

validate_pretty_results <- function(results_lin, endog_vars) {
  irf_array <- results_lin$irf_lin_mean

  # Select a test case
  test_impulse <- 2  # Example fixed index (adjustable)
  test_response <- 1
  test_horizon <- 6

  # Expected value from original array
  expected_value <- irf_array[test_response, test_horizon, test_impulse]

  # Run pretty_results and extract the full data frame
  irf_df <- pretty_results(results_lin, endog_vars)

  # Ensure impulse and response columns are character for comparison
  irf_df$impulse <- as.character(irf_df$impulse)
  irf_df$response <- as.character(irf_df$response)

  # Convert `endog_vars[test_impulse]` and `endog_vars[test_response]` to character for matching
  test_impulse_name <- as.character(endog_vars[test_impulse])
  test_response_name <- as.character(endog_vars[test_response])

  # Extract matching row
  matched_row <- subset(irf_df, 
                        impulse == test_impulse_name & 
                        response == test_response_name & 
                        horizon == test_horizon)

  # Print results
  print(paste("Expected IRF Value:", expected_value))

  # Ensure there's at least one match before printing
  if (nrow(matched_row) == 1) {
    print(paste("Extracted IRF Value:", matched_row$irf_value))

    if (round(expected_value, 6) == round(matched_row$irf_value, 6)) {
      print("✅ Validation passed!")
    } else {
      print("❌ Validation failed!")
    }
  } else {
    print("❌ No matching row found!")
  }
}
