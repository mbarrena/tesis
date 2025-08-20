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
#' @param fixed_lags An integer indicated a fixed number of lags. If specified, superseeds max_lags and lags_criterion.
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
run_lp_model <- function(
    df, 
    endog, 
    exog=NULL, 
    max_lags, 
    lags_criterion = 'AIC', 
    fixed_lags = NA, 
    newey_lags = NULL, 
    horizons = 10, 
    signif, 
    lags_exog = NULL, 
    trend = 0, 
    cumulative = FALSE, 
    threshold_var = NULL,
    upper_threshold = 0.80,
    lower_threshold = NULL,#0.20,
    discard_threshold = NULL#or 'upper' #or 'lower'   
  ) {
  
  data <- df
  dummy <- NULL

  if (!is_null(threshold_var)) {
    res <- prepare_threshold_dummy(df, threshold_var, upper_threshold, lower_threshold, discard_threshold)
    data <- res$df
    dummy <- res$dummy
  }

  # Map confidence levels to the corresponding values
  confint_map <- c("0.95" = 1.96, "0.68" = 1)

  # Check if the provided significance level is valid
  if (!as.character(signif) %in% names(confint_map)) {
    stop("Invalid significance level. Use 0.05 for 95% or 0.32 for 68%.")
  }

  if (!is.null(exog)) {
    exog_data <- data[, exog, drop = FALSE]
    if (is.null(lags_exog)) {
      if  (!is.na(fixed_lags)) {
        lags_exog=fixed_lags
      } else {
        lags_exog = max_lags
      }
    }
  } else {
    exog_data = NULL
  }

  if (!is.na(fixed_lags)) {
    max_lags = NA
    lags_criterion = NaN
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
      lags_criterion = lags_criterion,
      max_lags = max_lags,
      lags_endog_lin = fixed_lags,
      trend          = trend,  
      shock_type     = 1,  
      confint        = confint_map[as.character(signif)], 
      nw_lag         = newey_lags,
      nw_prewhite    = TRUE, 
      hor            = horizons,
      lags_exog = lags_exog
    )
  } else {
    switching <- dummy
    use_logistic <- TRUE
    if (!all(switching %in% c(0, 1))) {
      display_html("Dummy variable must have only 0 and 1 values")
      stop()
    }
    unique_count <- length(unique(switching))

    results_lin <- lp_nl(
      endog_data, 
      exog_data      = exog_data,
      lags_criterion = lags_criterion,
      max_lags = max_lags,
      lags_endog_lin = fixed_lags,
      lags_endog_nl = fixed_lags,
      trend          = trend,  
      shock_type     = 1,  
      confint        = confint_map[as.character(signif)], 
      nw_lag         = newey_lags,
      nw_prewhite    = TRUE, 
      hor            = horizons,
      lags_exog = lags_exog,
      switching = switching,
      use_logistic = FALSE,
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
    if (is_null(threshold_var)) {
      for (i in 1:dim(results_lin$irf_lin_mean)[3]) {
        # Cumsum mean IRFs
        results_lin$irf_lin_mean[,,i] <- t(apply(results_lin$irf_lin_mean[,,i], 1, cumsum))
        
        # Cumsum error bands to match cumulative mean
        results_lin$irf_lin_low[,,i] <- t(apply(results_lin$irf_lin_low[,,i], 1, cumsum))
        results_lin$irf_lin_up[,,i]  <- t(apply(results_lin$irf_lin_up[,,i], 1, cumsum))
      }
    } else {
      for (i in 1:dim(results_lin$irf_s1_mean)[3]) {
      # Cumsum mean IRFs
      results_lin$irf_s1_mean[,,i] <- t(apply(results_lin$irf_s1_mean[,,i], 1, cumsum))
      
      # Cumsum error bands to match cumulative mean
      results_lin$irf_s1_low[,,i] <- t(apply(results_lin$irf_s1_low[,,i], 1, cumsum))
      results_lin$irf_s1_up[,,i]  <- t(apply(results_lin$irf_s1_up[,,i], 1, cumsum))
    }
    for (i in 1:dim(results_lin$irf_s2_mean)[3]) {
      # Cumsum mean IRFs
      results_lin$irf_s2_mean[,,i] <- t(apply(results_lin$irf_s2_mean[,,i], 1, cumsum))
      
      # Cumsum error bands to match cumulative mean
      results_lin$irf_s2_low[,,i] <- t(apply(results_lin$irf_s2_low[,,i], 1, cumsum))
      results_lin$irf_s2_up[,,i]  <- t(apply(results_lin$irf_s2_up[,,i], 1, cumsum))
    }
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
    pretty_results_nl(results_lin = results_lin, endog_vars = endog, unique_count)
    # Process and print results with cumulative option
    plotlp_nl(results_lin, endog = endog, title_text = title_text, count_sections = unique_count)
  }

  return(results_lin)
}

prepare_threshold_dummy <- function(df, threshold_var, upper_threshold, lower_threshold, discard_threshold) {
  data <- df
  if (is.null(upper_threshold) && is.null(lower_threshold)) {
    display_html("Si se va a usar un umbral, especificar límite superior (con upper_threshold) o inferior (con lower_threshold)")
    stop()
  }
  if (! is.null(discard_threshold)) {
    # Caso thresholds descartando datos a partir de uno de ellos
    if (is.null(upper_threshold) || is.null(lower_threshold)) {
      display_html("Si se va a descartar datos en base a umbrales (discard_threshold = NULL), especificar tanto límite superior (con upper_threshold) como inferior (con lower_threshold)")
      stop()
    }
    upper_p <- quantile(data[[threshold_var]], upper_threshold, na.rm = TRUE)
    lower_p <- quantile(data[[threshold_var]], lower_threshold, na.rm = TRUE)

    df$above_upper_p <- ifelse(df[[threshold_var]] > upper_p, 1, 0)
    df$not_below_lower_p <- ifelse(df[[threshold_var]] < lower_p, 0, 1)
    if (discard_threshold == 'upper') {
      display_html("Se descartan los datos por arriba del umbral superior")
      data <- df[df$above_upper_p == 0, ]
      rownames(data) <- NULL
      dummy <- ifelse(data[[threshold_var]] < lower_p, 0, 1)
    } else if (discard_threshold == 'lower') {
      display_html("Se descartan los datos por debajo del umbral inferior")
      data <- df[df$not_below_lower_p == 1, ]
      rownames(data) <- NULL
      dummy <- ifelse(data[[threshold_var]] > upper_p, 1, 0)
    } else {
      display_html("discard_threshold debe ser 'upper' o 'lower'")
      stop()
    }
  } else {
    # Caso thresholds sin descartar datos
    if (! is.null(upper_threshold)) {
      upper_p <- quantile(data[[threshold_var]], upper_threshold, na.rm = TRUE)
      data$above_upper_p <- ifelse(data[[threshold_var]] > upper_p, 1, 0)
      dummy <- data$above_upper_p
      lower_threshold <- NULL
    }
    else if (! is.null(lower_threshold)) {
      lower_p <- quantile(data[[threshold_var]], lower_threshold, na.rm = TRUE)
      data$not_below_lower_p <- ifelse(data[[threshold_var]] < lower_p, 0, 1)
      dummy <- data$not_below_lower_p
      upper_threshold <- NULL
    }
  }
  
  plot_thresholds(df, threshold_var=threshold_var, upper_threshold=upper_threshold, lower_threshold=lower_threshold, discard_threshold=discard_threshold)
  return(list(df = data, dummy = dummy))
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

plotlp_nl <- function(results_lin, endog, title_text, count_sections,
                     plot_base_width = 5, plot_height = 4) {
  # Get the linear plots
  plot_width <- plot_base_width * count_sections
  linear_plots <- plot_nl(results_lin)
  n_sections <- length(linear_plots)
  
  # Create list to store plots
  all_plots <- list()
  
  # Set up the graphics device for notebook display
  # This needs to happen BEFORE creating plots
  options(repr.plot.width = plot_width, repr.plot.height = plot_height)
  
  # Loop through each variable pair
  for (i in seq_along(linear_plots[[1]])) {
    # Extract all sections' plots for this pair
    plot_list <- lapply(linear_plots, function(section) section[[i]])
    
    # Expand plot widths by modifying the theme
    plot_list <- lapply(plot_list, function(p) {
      p + theme(
        plot.margin = margin(5, 20, 5, 20),  # Add margins (top, right, bottom, left)
        aspect.ratio = NULL,  # Remove aspect ratio constraints
        panel.spacing = unit(2, "cm")  # Increase spacing between panels if faceted
      )
    })
    
    # Combine horizontally with dynamic width
    combined_plot <- cowplot::plot_grid(
      plotlist = plot_list,
      nrow = 1,
      labels = paste("Section", 1:n_sections),
      label_size = 10,
      rel_widths = rep(1, n_sections),  # Equal width for each section
      align = "h"  # Align horizontally
    )
    all_plots[[i]] <- combined_plot
  }
  
  # Print each plot
  for (i in seq_along(all_plots)) {
    if (!missing(title_text) && i == 1) {
      title_plot <- cowplot::ggdraw() +
        cowplot::draw_label(title_text, fontface = 'bold', size = 14)
      
      final_plot <- cowplot::plot_grid(
        title_plot,
        all_plots[[i]],
        ncol = 1,
        rel_heights = c(0.1, 1)
      )
    } else {
      final_plot <- all_plots[[i]]
    }
    
    # Force specific dimensions for display
    if (requireNamespace("gridExtra", quietly = TRUE)) {
      # gridExtra approach for better control
      grid::grid.newpage()
      g <- gridExtra::grid.arrange(final_plot, 
                                  widths = unit(plot_width, "in"), 
                                  heights = unit(plot_height, "in"))
      invisible(g)
    } else {
      # Regular print if gridExtra not available
      print(final_plot)
    }
  }
  
  invisible(all_plots)
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

pretty_results_nl <- function(results_lin, endog_vars, count_sections) {
  # Create a list to store results for each section
  all_sections_dfs <- list()
  
  # Process each section
  for (section in 1:count_sections) {
    # Extract IRF values, lower and upper bounds for this section
    irf_array <- results_lin[[paste0("irf_s", section, "_mean")]]
    lower_bound_array <- results_lin[[paste0("irf_s", section, "_low")]]
    upper_bound_array <- results_lin[[paste0("irf_s", section, "_up")]]
    
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
      section = paste("Section", section),
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
    
    # Map impulse and response numbers to names
    irf_df$impulse <- endog_vars[irf_df$impulse]
    irf_df$response <- endog_vars[irf_df$response]
    
    # Add to the list of all sections
    all_sections_dfs[[section]] <- irf_df
  }
  
  # Combine all sections into one data frame
  combined_df <- do.call(rbind, all_sections_dfs)
  
  # Group by impulse-response pair (regardless of section)
  var_pairs <- unique(combined_df[, c("impulse", "response")])
  
  # For each variable pair, create a multi-index dataframe with all sections
  for (i in 1:nrow(var_pairs)) {
    impulse_var <- var_pairs$impulse[i]
    response_var <- var_pairs$response[i]
    
    # Filter data for this variable pair
    pair_data <- combined_df[combined_df$impulse == impulse_var & 
                             combined_df$response == response_var, ]
    
    # Create a wider format with sections as columns
    # First pivot for irf_value
    wide_irf <- reshape(pair_data, 
                        idvar = "horizon",
                        timevar = "section", 
                        direction = "wide",
                        v.names = c("irf_value", "lower_bound", "upper_bound"))
    
    # Sort by horizon
    wide_irf <- wide_irf[order(wide_irf$horizon), ]
    
    # Display the header
    header_html <- paste("<h3>Impulse:", impulse_var, "-> Response:", response_var, "</h3>")
    display_html(header_html)
    
    # Rename columns to be more readable
    names(wide_irf) <- gsub("irf_value.Section ", "IRF S", names(wide_irf))
    names(wide_irf) <- gsub("lower_bound.Section ", "Lower S", names(wide_irf))
    names(wide_irf) <- gsub("upper_bound.Section ", "Upper S", names(wide_irf))
    
    # Remove redundant columns
    wide_irf <- wide_irf[, !names(wide_irf) %in% c("impulse", "response")]
    
    # Display the multi-section dataframe
    display(wide_irf)
  }
  
  return(combined_df)
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
    lags_endog_nl = if (!is.null(specs$lags_endog_nl)) specs$lags_endog_nl else NA,
    lags_criterion = specs$lags_criterion,
    max_lags = specs$max_lags,
    lags_exog = if (!is.null(specs$lags_exog)) specs$lags_exog else NA,
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

plot_thresholds <- function(df, threshold_var, upper_threshold=NULL, lower_threshold=NULL, discard_threshold=NULL) {
  par(bg = "white") # mar increases bottom margin for taller labels

  # Create a line plot of ipc with x-axis as the dataframe index
  plot(df$ipc, type = "l", col = "black", lwd = 2,
      xlab = "Time (Year + Quarter)", ylab = "IPC", main = "IPC with thresholds", xaxt = "n")

  # Add custom x-axis labels
  ticks <- seq(4, nrow(df), by = 4)
  axis(1, at = ticks, labels = paste(df$año[ticks], df$trimestre[ticks]), las = 2, cex.axis = 0.7)

  if (is.null(discard_threshold)) {
    color_upper = "blue"
    color_lower = "blue"
  } else {
    remove_lower = discard_threshold == 'lower'
    color_upper = ifelse(remove_lower, "green", "red")
    color_lower = ifelse(remove_lower, "red", "green")
  }

  legend = c('IPC')
  colors = c("black")
  # Calculate the 80th percentile and plot
  if (!is.null(upper_threshold)) {
    upper_p <- quantile(df[[threshold_var]], upper_threshold, na.rm = TRUE)
    abline(h = upper_p, col = color_upper, lwd = 2, lty = 2)
    legend <- append(legend, paste0(upper_threshold*100,"th Percentile (", round(upper_p, 4), ")"))
    colors <- append(colors, color_upper)
    
    # Add shading for upper threshold if it's the discard threshold
    if (!is.null(discard_threshold) && discard_threshold == 'upper') {
      # Get plot limits
      y_lim <- par("usr")[3:4]
      x_lim <- par("usr")[1:2]
      
      # Create polygon coordinates for the area above the threshold
      x_coords <- c(x_lim[1], x_lim[2], x_lim[2], x_lim[1])
      y_coords <- c(upper_p, upper_p, y_lim[2], y_lim[2])
      
      # Add semi-transparent red polygon
      polygon(x_coords, y_coords, col = rgb(1, 0, 0, 0.2), border = NA)
      
      # Add to legend
      legend <- append(legend, "Discarded Region")
      colors <- append(colors, rgb(1, 0, 0, 0.2))
    }
  }
  
  if (!is.null(lower_threshold)) {
    lower_p <- quantile(df[[threshold_var]], lower_threshold, na.rm = TRUE)
    abline(h = lower_p, col = color_lower, lwd = 2, lty = 2)
    legend <- append(legend, paste0(lower_threshold*100,"th Percentile (", round(lower_p, 4), ")"))
    colors <- append(colors, color_lower)
    
    # Add shading for lower threshold if it's the discard threshold
    if (!is.null(discard_threshold) && discard_threshold == 'lower') {
      # Get plot limits
      y_lim <- par("usr")[3:4]
      x_lim <- par("usr")[1:2]
      
      # Create polygon coordinates for the area below the threshold
      x_coords <- c(x_lim[1], x_lim[2], x_lim[2], x_lim[1])
      y_coords <- c(y_lim[1], y_lim[1], lower_p, lower_p)
      
      # Add semi-transparent red polygon
      polygon(x_coords, y_coords, col = rgb(1, 0, 0, 0.2), border = NA)
      
      # Add to legend
      legend <- append(legend, "Discarded Region")
      colors <- append(colors, rgb(1, 0, 0, 0.2))
    }
  }

  # Add a legend to the plot
  legend("topright", legend = legend,
        col = colors, lty = c(1, rep(2, length(colors)-2), NA), lwd = c(2, rep(2, length(colors)-2), NA), 
        fill = c(NA, rep(NA, length(colors)-2), rgb(1, 0, 0, 0.2)))
}