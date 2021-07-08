#' Simulate covariates
#'
#'
#' @export
#'
pdSimuX <- function(n_pat, muCov, sdCov, corCov, cov.breaks = NULL,
                    fmla_x = NULL, coeff_x = 1, dx = NULL, seed = NULL, ...) {

    fmla <-  fmla_x;
    stopifnot(n_pat > 0 & is.null(dx));
    stopifnot(inherits(fmla,"formula") | is.null(fmla));

    if (!is.null(seed))
        old.seed <- set.seed(seed);

    ## simulate covariates
    if (is.null(dx)) {
        cov_x <- rmvnorm(n_pat,
                         mean  = muCov,
                         sigma = get.covmat(sdCov, corCov));

        colnames(cov_x) <- paste("V", 1:ncol(cov_x), sep = "");
        cov_x           <- data.frame(cov_x);
        cov_x           <- get.cov.cat(cov_x, cov.breaks);
    }

    ## xbeta
    if (is.null(fmla)) {
        fmla <- formula(paste("~ -1 + ",
                              paste(colnames(cov_x), collapse = "+")));
    }

    d.matrix <- model.matrix(fmla, cov_x);
    xbeta    <- get.xbeta(d.matrix, coeff_x);

    ## return
    if (!is.null(seed))
        set.seed(old.seed)

    list(xbeta = xbeta,
         dx    = cov_x);
}

#' Simulate enrollment
#'
#'
#' @export
#'
pdSimuTime <- function(n_pat, fmla_t = NULL, coeff_t = 1,
                       dtime = NULL,
                       t_start = 0, t_end = 3, seed = NULL, ...) {

    ## check par
    fmla <-  fmla_t;
    stopifnot(inherits(fmla, "formula") | is.null(fmla));

    if (!is.null(seed))
        old_seed <- set.seed(seed);


    ## simulate time
    if (is.null(dtime)) {
        ## rexp(n = n_pat, rate = exp_rate)
        dtime <- data.frame(time = runif(n_pat, min = t_start, max = t_end))
    }

    ## tgamma
    if (is.null(fmla)) {
        fmla <- formula(paste("~ -1 + time"));
    }

    d_matrix <- model.matrix(fmla, dtime);
    tgamma   <- get.xbeta(d_matrix, coeff_t);

    if (!is.null(seed))
        set.seed(old_seed)

    list(tgamma = tgamma,
         dtime  = dtime);
}

#' Simulate enrollment
#'
#'
#' @export
#'
pdSimuPts <- function(n_pat, mu_0 = 0, ...) {

    ## time
    simu_t <- pdSimuTime(n_pat, ...);

    ## covariates
    simu_x <- pdSimuX(n_pat, ...);

    ## outcome
    mu_all <- mu_0 + simu_t$tgamma + simu_x$xbeta;
    emu    <- expit(mu_all);
    y      <- rbinom(n_pat, 1, emu);

    ## return
    cbind(data.frame(y = y, xbeta = simu_x$xbeta,
                     tgamma = simu_t$tgamma, mu = emu),
          simu_x$dx, simu_t$dtime);
}

#' Create time interval
#'
#'
#' Map enrollment time to [0, start of current study] and [start of current
#' study, end of current study]. Then, create intervals for [0, start of current
#' study].
#'
#' @param vec_t patient enrollment time
#' @param trange map enrollment time to certain range
#' @param n_interval number of intervals for historical enrollment
#' @param t_start time corresponding to the start of the study
#' @param t_hist time corresponding to the end of the historical study and start
#'               of the current randomized segment
#' @param t_end time corresponding to the end of the study
#'
#' @export
#'
pdGetInt  <- function(vec_t, n_interval = 4, t_start = 0,
                      t_hist = 2, t_end = 3, remap = FALSE, ...) {

    map_t <- vec_t
    if (remap) {
        ## map time
        max_t <- max(vec_t) + 0.0001
        min_t <- min(vec_t) - 0.0001
        lt    <- t_end - t_start
        map_t <- t_start + lt * (vec_t - min_t) / (max_t - min_t)
    }

    ## intervals
    lint  <- (t_hist - t_start) / n_interval
    vec_i <- ceiling(map_t / lint)
    vec_i[which(map_t > t_hist)] <- n_interval + 1

    ## interval time
    itime <- c(seq(lint/2, t_hist - lint/2, lint),
               t_hist + (t_end - t_hist)/2)

    ## label
    sinterval <- c(paste("Interval", 1:n_interval),  "Concurrent")

    ## return
    data.frame(orgt     = vec_t, ## original time
               mapt     = map_t,
               interval = vec_i,
               itime    = itime[vec_i],
               ilabel   = sinterval[vec_i])
}
