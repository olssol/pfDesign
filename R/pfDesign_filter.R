## -------------------------------------------------------------------------
##
##               FUNCTIONS FOR PARTICLE FILTERING ALGORITHMS
##
## -------------------------------------------------------------------------

#' Particle filter all
#'
#' @param vecy vector of binary outcomes
#' @param mix_ab beta shape parameters (a, b) for a non-informative prior
#' @param epsillon mixture proportion from the prior samples. ignored if
#'     prior_smps is NULL
#' @param f_tran transition function
#' @param f_post posterior sampling when the prior is pi_0, usually for the
#'     first interval
#' @param f_pi0 pi_0 function
#' @param f_filter particle filter function
#' @param update_eps whether keep epsilon as constant (design parameter)
#' @param ... parameters for the transition function
#'
#' @export
#'
pdFilter <- function(vec_y, vec_interval, nsmps = 5000, epsilons = 1,
                     f_tran = pdTran, f_post = pdBinPost, f_pi0 = pdBinPi0,
                     f_filter = pdFilterSingle, update_eps = FALSE,
                     ...) {

    f_d <- function(p, l, i, eps) {
        rst  <- rbind(data.frame(theta = as.numeric(p),
                                 type  = "Posterior"),
                      data.frame(theta = as.numeric(l),
                                 type  = "Mixed"))
        rst$interval <- i
        rst$epsilon  <- eps
        rst
    }

    ## check intervals
    intvs  <- unique(vec_interval)
    n_intv <- length(intvs)

    stopifnot(1      == min(intvs) &
              n_intv == max(intvs));

    ## fill epsilons
    if (length(epsilons) < n_intv) {
        epsilons <- rep(epsilons, n_intv);
    }

    ## samples after the first interval
    y1           <- vec_y[which(1 == vec_interval)]
    post_smps    <- f_post(y1, nsmps = nsmps, ...)
    last_weights <- rep(1 / nsmps, nsmps)

    ## transition
    last_eps  <- epsilons[1]
    last_smps <- f_tran(post_smps, epsilon = last_eps, f_pi0, ...)
    all_rst   <- f_d(post_smps, last_smps[, 1], 1, eps = last_eps)

    ## rest intervals
    if (n_intv > 1) {
        for (i in 2:n_intv) {
            cury <- vec_y[which(i == vec_interval)]

            ## epsilon
            if (!update_eps) {
                last_eps <- epsilons[i]
            } else {
                ## to be added
                ## last_epsilon <- mean(posts[, 2])
            }

            ## posterior sampling
            posts     <- f_filter(cury,
                                  last_smps    = last_smps,
                                  last_weights = last_weights,
                                  nsmps        = nsmps,
                                  ...)
            post_smps <- posts[, 1, drop = FALSE]

            ## transition
            last_smps    <- f_tran(post_smps, epsilon = last_eps, f_pi0, ...)
            last_weights <- posts[, 3]


            ## append results
            all_rst <- rbind(all_rst,
                             f_d(post_smps, last_smps[, 1], i, last_eps))
        }
    }

    ## return
    all_rst
}


#' Particle filter for a given interval
#'
#'
#' @param y current Y
#' @param last_smps col 1: prior samples of theta; col 2: mixture indicator
#' @param last_weights normalized weights corresponding to last_smps. If NULL,
#'     all weights are 1/N.
#' @param f_ll likelihood function
#'
#' @return col 1: theta; col 2: mixture indicator; col 3: weights
#' @export
#'
pdFilterSingle <- function(y, last_smps, last_weights = NULL,
                           f_ll = pdBinPmf, thresh_ess = NULL, ...) {

    ns <- nrow(last_smps)

    ## weights
    if (is.null(last_weights))
        last_weights <- 1 / n

    weights <- f_ll(last_smps[, 1], y)
    weights <- weights * last_weights
    cur_w   <- weights / sum(weights)

    ##  samples
    cur_smps <- last_smps

    if (!is.null(thresh_ess)) {
        ess      <- 1 / sum(cur_w^2)
        do_resmp <- ess < thresh_ess
    } else {
        do_resmp <- TRUE
    }

    if (do_resmp) {
        next_inx <- sample(1:ns, replace = TRUE, prob = cur_w)
        cur_smps <- cur_smps[next_inx,]
        cur_w    <- 1 / ns
    }

    ## return
    cbind(cur_smps, cur_w)
}


#' HMM transition function
#'
#' @inheritParams parameters
#'
#'
#' @return column 1: theta; column 2: mixture indicator
#'
#' @export
#'
pdTran <- function(smps, epsilon = 1, f_pi0 = pdBinPi0, ...) {
    n    <- nrow(smps)
    u    <- rbinom(n, 1, epsilon)
    frm0 <- which(0 == u)
    n0   <- length(frm0)

    if (n0 > 0) {
        smp0       <- f_pi0(n0, ...)
        smps[frm0] <- smp0
    }

    ## return
    cbind(smps, u)
}


#' Binomial likelihood function
#'
#'
#'
#' @export
#'
pdBinPmf <- function(theta, y) {
    theta^sum(y) * (1 - theta)^sum(1 - y)
}



#' Noninformative Pi_0(theta)
#'
#'
#'
#' @inheritParams parameters
#'
#'
#' @export
#'
pdBinPi0 <- function(nsmps, mix_ab = c(0.25, 0.25), ...) {
    rbeta(nsmps,
          mix_ab[1],
          mix_ab[2])
}


#' Posterior Beta
#'
#' @inheritParams parameters
#'
#' @export
#'
pdBinPost <- function(y, nsmps = 5000, mix_ab = c(0.25, 0.25), ...) {
    rst <- rbeta(nsmps,
                 sum(y)     + mix_ab[1],
                 sum(1 - y) + mix_ab[2])
    cbind(rst)
}
