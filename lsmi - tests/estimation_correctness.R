# checking if estimated MI values correspond to real ones

library(easypackages)
packages("mvtnorm", "dplyr", "ggplot2")
source("./lsmi - method/LSMI_extra.R")

set.seed(451)
# simulating big samples of multivariate normal with
# various correlation coefficients to compare against true MI values

# taken from Wikipedia
mutual_info_2dN <- function(rho) -0.5 * log2(1 - rho^2)
# generates a n * 2 matrix of bivariate normal observations
generate_2dN_sample <- function(rho, n) rmvnorm(n, sigma = matrix(c(1, rho, rho, 1), nrow=2))

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

compare_sample_size <- function(rhos, ns, lsmi_pars) {
    # in: correlation coefficients, vector of sample sizes, number of repeats, lsmi estimator parameters
    # out: df with cols (mutual_info, lsmi_estimate, sample_x, sample_y, sample_size)
    
}

compare_n_base_functions <- function(rhos, n, lsmi_pars) {

}