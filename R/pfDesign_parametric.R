#' Estimate response rates using parametric logistic model
#'
#' @param dta  data frame
#' @param fml  logistic model formula
#' @param n_bs number of bootstraps
#'
#' @return mean and bootstrap variance of the current study response rates
#'
#' @export
#'
pdLogist <- function(dta, fml, n_bs = 500) {

    f_est <- function(d, intervals) {
        r_glm <- glm(fml, data = d, family = "binomial")
        r_prd <- predict(r_glm, type = "response")

        rst <- NULL
        for (i in intervals) {
            cinx <- which(dta$interval == i)
            if (0 == length(cinx)) {
                rst <- c(rst, NA)
            } else {
                rst <- c(rst, mean(r_prd[cinx]))
            }
        }
        rst
    }

    ## unique intervals
    intervals <- sort(unique(dta$interval))

    ## original data
    m <- f_est(dta, intervals)

    ## bootstrap
    vec_bs <- NULL
    for (i in 1:n_bs) {
        cur_b  <- sample(seq_len(nrow(dta)), replace = TRUE)
        cur_d  <- dta[cur_b, ]
        vec_bs <- cbind(vec_bs,
                        f_est(cur_d, intervals))
    }

    cbind(interval = intervals,
          mean = m,
          variance = apply(vec_bs, 1, var, na.rm = TRUE))
}
