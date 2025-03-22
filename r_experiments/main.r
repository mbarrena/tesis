
run_lp_model <- function(data, endog, max_lags, newey_lags = NULL, horizons = 10, signif) {
  # Map confidence levels to the corresponding values
  confint_map <- c("0.05" = 1.96, "0.32" = 1)
  
  # Check if the provided significance level is valid
  if (!as.character(signif) %in% names(confint_map)) {
    stop("Invalid significance level. Use 0.05 for 95% or 0.32 for 68%.")
  }
  
  # Select endogenous variables
  endog_data <- data[, endog, drop = FALSE]
  
  # Convert to numeric (ensure proper model fitting)
  endog_data <- data.frame(lapply(endog_data, as.numeric))
  
  # Run the local projections model
  results_lin <- lp_lin(
    endog_data, 
    lags_endog_lin = max_lags,  
    trend          = 0,  
    shock_type     = 1,  
    confint        = confint_map[as.character(signif)], 
    nw_lag         = newey_lags,
    hor            = horizons
  )
  
  # Generate and return plots
  linear_plots <- plot_lin(results_lin)
  
  # Show all plots
  lin_plots_all <- sapply(linear_plots, ggplotGrob)
  marrangeGrob(lin_plots_all, nrow = length(endog), ncol = length(endog), top = NULL)
}

# Example of how to call the function
# Assuming `data` is already read from the Excel file
run_lp_model(data, c("impp_usa", "E", "ipc_ajust", "pbird"), max_lags = 3, newey_lags = 4, signif = 0.05)
