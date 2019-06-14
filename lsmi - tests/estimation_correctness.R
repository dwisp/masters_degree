# checking if estimated MI values correspond to real ones

library(easypackages)
packages("mvtnorm", "dplyr", "ggplot2")
source("./lsmi - method/LSMI_extra.R")

set.seed(451)
# simulating big samples of multivariate normal with
# various correlation coefficients to compare against true MI values

# taken from Wikipedia
mutual_info_2dN <- function(rho) -0.5 * log2(1 - rho^2)

# custom example
mutual_info_Nexp <- function(lambda, sigma) {
    mi <- -0.5 - log(sqrt(2 * pi) * sigma)
    f_y <- function(y, s = sigma) pnorm((1 - y) / s) -  pnorm((- y / s))
    y_entr <- function(y) f_y(y) * log(f_y(y))
    integral <- integrate(y_entr, -Inf, Inf)
    return(mi - integral[["value"]])
}

# generates a n * 2 matrix of bivariate normal observations
generate_2dN_sample <- function(rho, n, seed=NULL) {
    if(!is.null(seed)) set.seed(seed)
    return(rmvnorm(n, sigma = matrix(c(1, rho, rho, 1), nrow=2)))
}

# custom more complicated example
generate_2dNexp_sample <- function(lambda, sigma, n) {
    x_s <- rexp(n, lambda)
    # p(y|x) is defined through exponentiated x
    y_s <- rnorm(n, exp(-lambda * x_s), sigma)
    return(cbind(x_s, y_s))
}

generate_sample_df <- function(rhos=seq(0, 0.9, 0.1), n=1000) {
    # list of correlation matrices
    # sigmas <- lapply(rhos, function(rho) matrix(1, rho, rho, 1))
    # generating samples of given size
    samples <- lapply(rhos, function(rho) generate_2dN_sample(rho, n))
    samples_x <- lapply(samples, function(smpl) smpl[, 1])
    samples_y <- lapply(samples, function(smpl) smpl[, 2])

    true_mi_values <- sapply(rhos, mutual_info_2dN)

    df <- tibble(
        "mutual_info" = true_mi_values,
        "lsmi_estimate" = 0,
        "sample_x" = samples_x,
        "sample_y" = samples_y,
        # "lsmi_params" = list()
    )
    return(df)
}

estimate_lsmi <- function(df, lsmi_pars=list("method.nbfuns" = "uniform")) {
    df[["lsmi_estimate"]] <- mapply(
        lsmi.extra,
        df[["sample_x"]],
        df[["sample_y"]],
        MoreArgs = lsmi_pars,
        SIMPLIFY = TRUE
    )
    return(df)
}

compare_sampling_methods <- function(rhos, n, cnt, lsmi_pars=list()) {
    # in: correlation coefficients, sample size, number of repeats, lsmi estimator parameters
    # out: df with cols (mutual_info, lsmi_estimate, sample_x, sample_y, base_funs)
    df_uniform <-
        generate_sample_df(rhos, n) %>%
        mutate(base_funs="uniform")
    df_uniform <- df_uniform[rep(1:nrow(df_uniform), each=cnt), ]

    df_nonpaired <-
        generate_sample_df(rhos, n) %>%
        mutate(base_funs="non-paired")
    df_nonpaired <- df_nonpaired[rep(1:nrow(df_nonpaired), each=cnt), ]

    df_uniform <- estimate_lsmi(df_uniform, c(lsmi_pars, "method.nbfuns" = "uniform"))
    df_nonpaired <- estimate_lsmi(df_nonpaired, c(lsmi_pars, "method.nbfuns" = "non-paired"))

    return(bind_rows(df_uniform, df_nonpaired))
}

compare_sample_size <- function(rhos, ns, lsmi_pars=list("method.nbfuns" = "uniform")) {
    # in: correlation coefficients, vector of sample sizes, lsmi estimator parameters
    # out: df with cols (mutual_info, lsmi_estimate, sample_x, sample_y, sample_size)
    dfs <- list()
    for (n in ns) {
        dfs[[as.character(n)]] <- generate_sample_df(rhos, n)
    }
    df_samples <- bind_rows(dfs)
    # sapply?
    df_samples %<>% mutate(sample_size = sapply(sample_x, length))
    df_samples <- estimate_lsmi(df_samples, lsmi_pars)
    return(df_samples)
}

compare_n_base_functions <- function(rhos, n, nbfuns, lsmi_pars=list("method.nbfuns" = "uniform")) {
    # in: correlation coefficients, sample sizes, vector of number of base functions, lsmi estimator parameters
    # out: df with cols (mutual_info, lsmi_estimate, sample_x, sample_y, sample_size)
    # stopifnot(!"nbfuns" %in% names(lsmi_pars), "Number of base functions is provided via nbfuns argument")
    # a single base data frame
    df <- generate_sample_df(rhos, n)
    results <- c()
    for (nbf in nbfuns) {
        tmp_df <- estimate_lsmi(df, c(lsmi_pars, "nbfuns" = nbf))
        results <- c(results, tmp_df[["lsmi_estimate"]])
    }
    df <- df[rep(1:nrow(df), length(nbfuns)), ]
    df %<>%
        mutate(
            lsmi_estimate = results,
            num_base_funs = rep(nbfuns, each = length(rhos))
        )
    return(df)
}
