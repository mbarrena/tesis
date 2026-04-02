fitSTVAR(
    data,
    p,
    M,
    weight_function = c(
        "relative_dens", "logistic", "mlogit", "exponential", "threshold",
        "exogenous"
    ), # Aca logistic para que sea VLSTAR
    weightfun_pars = NULL, # Si es logistic, vector c(el numero de la variable de transición, el numero de lag)
    cond_dist = c("Gaussian", "Student", "ind_Student", "ind_skewed_t"), # Nico recomienda probar Gaussiana o Student
    parametrization = c("intercept", "mean"),
    AR_constraints = NULL,
    mean_constraints = NULL,
    weight_constraints = NULL,
    estim_method,
    penalized,
    penalty_params = c(0.05, 0.2),
    allow_unstab,
    min_obs_coef = 3,
    sparse_grid = FALSE,
    h = 0.001,
    nrounds,
    ncores = 2,
    maxit = 2000,
    seeds = NULL,
    print_res = TRUE,
    use_parallel = TRUE,
    calc_std_errors = TRUE,
    ...
)
