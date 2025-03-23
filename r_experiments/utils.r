library(lpirfs)
library(gridExtra)
library(ggplot2)
library(readxl)
library(httr)

library(reshape2)
library(purrr) 

run_lp_model <- function(data, endog, exog, max_lags, newey_lags = NULL, horizons = 10, signif, lags_exog = NULL, cumulative = FALSE) {
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
  
  # Run the local projections model
  results_lin <- lp_lin(
    endog_data, 
    exog_data = exog_data,
    lags_endog_lin = max_lags,  
    trend          = 0,  
    shock_type     = 1,  
    confint        = confint_map[as.character(signif)], 
    nw_lag         = newey_lags,
    hor            = horizons,
    lags_exog = lags_exog
  )

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

  plotlp(results_lin, endog = endog, title_text = title_text)
  
  # Process and print results with cumulative option
  pretty_results(results_lin = results_lin, endog_vars = endog)

  return(results_lin)
}

plotlp <- function(results_lin, endog, title_text) {
  # Generate plots
  linear_plots <- plot_lin(results_lin)

  print(title_text)
  
  # Show all plots
  lin_plots_all <- sapply(linear_plots, ggplotGrob)

  final_plot <- marrangeGrob(lin_plots_all, nrow = length(endog), ncol = length(endog), 
                             top = grid::textGrob(title_text, gp = grid::gpar(fontsize = 14, fontface = "bold")))

  print(final_plot)
}

pretty_results <- function(results_lin, endog_vars) {
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
  lapply(irf_list, print)
  return(irf_named_df)
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
