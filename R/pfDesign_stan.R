#' Call STAN models
#'
#'
#' @param chains STAN parameter. Number of Markov chainsm
#' @param iter STAN parameter. Number of iterations
#' @param warmup STAN parameter. Number of burnin.
#' @param control STAN parameter. See \code{rstan::stan} for details.
#' @param ... other options to call STAN sampling such as \code{thin},
#'     \code{algorithm}. See \code{rstan::sampling} for details.#'
#'
#'
#' @export
#'
pdSTAN <- function(lst.data, stan.mdl = "mixprior",
                   chains = 4, iter = 2000, warmup = 1000,
                   control = list(adapt_delta=0.95), cores = 4, ...) {

    stan.rst <- rstan::sampling(stanmodels[[stan.mdl]],
                                data    = lst.data,
                                chains  = chains,
                                iter    = iter,
                                warmup  = warmup,
                                cores   = cores,
                                control = control,
                                ...);

    stan.rst;
}
