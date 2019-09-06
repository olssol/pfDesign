#' Sample with a mixture prior for binary outcomes
#'
#' @param vecy vector of binary outcomes
#' @param iter number of thetas to draw when there is no prior samples.
#'     otherwise, draw the same number of thetas as in the prior samples
#' @param mix_ab beta shape parameters (a, b) for a non-informative prior
#' @param epsillon mixture proportion from the prior samples. ignored if
#'     prior_smps is NULL
#'
#' @param ... optional parameters for \code{\link{pdStan}}
#' @return samples of \eqn{theta}
#'
#' @export
#'
pdSampleSingle <- function(vec_y, prior_smps = NULL, epsilon = 1, mix_ab = c(0.25, 0.25),
                           nsmps = 5000, chains = 4, warmup = 1000, ...) {
    Y    = sum(vec_y);
    N    = length(vec_y);
    SMP  = prior_smps;

    ## no prior samples then take non-informative priors
    if (is.null(SMP)) {
        rst <- rbeta(nsmps, Y + mix_ab[1], N - Y + mix_ab[2]);
        return(rst);
    }

    ## draw samples using stan
    mix_ind = rbinom(length(SMP), 1, 1-epsilon);
    n_ninf  = sum(mix_ind);

    if (n_ninf > 0) {
        SMP[which(0 == mix_ind)] = rbeta(n_ninf, mix_ab[1], mix_ab[2]);
    }

    ## sort samples
    SMP  = c(-Inf, sort(SMP), Inf);
    NSMP = length(SMP);

    ## call stan
    post_theta = pdSTAN(list(Y = Y, N = N, NSMP = NSMP, SMP = SMP),
                        stan.mdl = "mixprior",
                        iter   = (nsmps/chains) + warmup,
                        warmup = warmup,
                        chains = chains,
                        ...);

    post_theta
}


#' Draw all samples
#'
#'
#' @return samples of \eqn{theta}
#'
#' @export
#'
pdSample <- function(vec_y, vec_interval, epsilons = 1, ...) {

    ## check intervals
    intvs = unique(vec_interval);

    ## number of intervals
    n_intv    = length(intvs);

    stopifnot(1      == min(intvs) &
              n_intv == max(intvs));

    ## fill epsilons
    if (length(epsilons) < (n_intv - 1)) {
        epsilons = rep(epsilons, n_intv - 1);
    }
    ## fix epsilon NA
    epsilons = c(NA, epsilons);

    ## posterior samples
    all_rst   = NULL;
    last_smps = NULL;
    for (i in 1:n_intv) {
        cury      = vec_y[which(i == vec_interval)];
        currst    = pdSampleSingle(cury, prior_smps = last_smps, epsilon = epsilons[i], ...);
        last_smps = extract(currst, par = "theta")$theta;
        all_rst   = rbind(all_rst,
                          data.frame(theta = last_smps,
                                     interval = i));
    }

    ## return
    all_rst
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
    vec_t = t_start + lt*(vec_t - min_t)/(max_t - min_t);
    vec_i = ceiling(vec_t / lint);

    data.frame(mapt     = vec_t,
               interval = vec_i);

}
