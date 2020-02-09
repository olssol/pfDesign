#' Estimate response rates using parametric logistic model
#'
#' @param dta  data frame
#' @param fml  logistic model formula
#' @param n_bs number of bootstraps
#' @param inx  current study indices
#'
#' @return mean and bootstrap variance of the current study response rates
#'
#' @export
#'
pdLogist <- function(dta, fml, inx, n_bs = 500) {

    f_est <- function(d) {
        r_glm <- glm(fml, data = d, family = "binomial")
        r_prd <- predict(r_glm)[inx]
        mean(r_prd)
    }

    ## original data
    m  <- f_est(dta)

    ## bootstrap
    vec_bs <- NULL
    for (i in 1:n_bs) {
        cur_b  <- sample(seq_len(nrow(dta)), replace = TRUE)
        cur_d  <- dta[cur_b, ]
        vec_bs <- c(vec_bs, f_est(cur_d))
    }

    c(mean = m, variance = var(vec_bs))
}
