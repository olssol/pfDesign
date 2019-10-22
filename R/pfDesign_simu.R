#' Simulate covariates
#'
#'
#' @export
#'
pdSimuX <- function(nPat, muCov, sdCov, corCov, cov.breaks = NULL,
                    fmla_x = NULL, coeff_x = 1, dx = NULL, seed = NULL, ...) {

    fmla = fmla_x;
    stopifnot(nPat > 0 & is.null(dx));
    stopifnot(inherits(fmla,"formula") | is.null(fmla));

    if (!is.null(seed))
        old.seed <- set.seed(seed);

    ## simulate covariates
    if (is.null(dx)) {
        cov_x <- rmvnorm(nPat,
                         mean  = muCov,
                         sigma = get.covmat(sdCov, corCov));

        colnames(cov_x) <- paste("V", 1:ncol(cov_x), sep="");
        cov_x           <- data.frame(cov_x);
        cov_x           <- get.cov.cat(cov_x, cov.breaks);
    }

    ## xbeta
    if (is.null(fmla)) {
        fmla <- formula(paste("~ -1 + ",
                              paste(colnames(cov_x), collapse = "+")));
    }

    d.matrix = model.matrix(fmla, cov_x);
    xbeta    = get.xbeta(d.matrix, coeff_x);

    ## return
    if (!is.null(old.seed))
        set.seed(old.seed)

    list(xbeta = xbeta,
         dx    = cov_x);
}

#' Simulate enrollment
#'
#'
#' @export
#'
pdSimuTime <- function(nPat, exp_rate = 1, fmla_t = NULL, coeff_t = 1, dtime = NULL, seed = NULL, ...) {

    ## check par
    fmla = fmla_t;
    stopifnot(inherits(fmla,"formula") | is.null(fmla));

    if (!is.null(seed))
        old.seed <- set.seed(seed);


    ## simulate time
    if (is.null(dtime)) {
        ## rexp(n = nPat, rate = exp_rate)
        dtime = data.frame(time = runif(nPat));
    }

    ## tgamma
    if (is.null(fmla)) {
        fmla <- formula(paste("~ -1 + time"));
    }

    d.matrix = model.matrix(fmla, dtime);
    tgamma   = get.xbeta(d.matrix, coeff_t);

    if (!is.null(seed))
        set.seed(old.seed)

    list(tgamma = tgamma,
         dtime  = dtime);
}

#' Simulate enrollment
#'
#'
#' @export
#'
pdSimuPts <- function(nPat, mu_0 = 0, ...) {

    ## time
    simu_t = pdSimuTime(nPat, ...);

    ## covariates
    simu_x = pdSimuX(nPat, ...);

    ## outcome
    mu_all = mu_0 + simu_t$tgamma + simu_x$xbeta;
    emu    = expit(mu_all);
    y      = rbinom(nPat, 1, emu);

    ## return
    cbind(data.frame(y = y, xbeta = simu_x$xbeta,
                     tgamma = simu_t$tgamma, mu = emu),
          simu_x$dx, simu_t$dtime);
}

#' Create interval
#'
#' @param trange map enrollment time to certain range
#'
#' @export
#'
pdGetInt  <- function(vec_t, n_interval = 4, t_start = 0, t_end = 2) {
    max_t = max(vec_t) + 0.0001;
    min_t = min(vec_t) - 0.0001;

    lt    = t_end - t_start;
    lint  = lt / n_interval;
    map_t = t_start + lt*(vec_t - min_t)/(max_t - min_t);
    vec_i = ceiling(map_t / lint);

    data.frame(mapt     = map_t,
               interval = vec_i);

}

