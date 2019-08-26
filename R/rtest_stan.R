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
rtSTAN <- function(lst.data, stan.mdl = "example",
                   chains = 4, iter = 2000, warmup = 1000,
                   control = list(adapt_delta=0.95), ...) {

    stan.rst <- rstan::sampling(stanmodels[[stan.mdl]],
                                data    = lst.data,
                                chains  = chains,
                                iter    = iter,
                                warmup  = warmup,
                                control = control,
                                ...);

    stan.rst;
}
