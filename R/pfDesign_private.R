##----------------------------------------------------------------------------
##                GENERAL FUNCTIONS
##----------------------------------------------------------------------------

make.global <- function(alist, dest.env='.GlobalEnv') {
    for (i in 1:length(alist)) {
        assign(names(alist[i]), alist[[i]], dest.env );
    }
}


expit <- function(x) {
    ex <- exp(x);
    ex/(1+ex);
}

get.xbeta <- function(covX, regCoeff) {
    if (length(regCoeff) > 0 &
        length(regCoeff) != ncol(covX))
        warning("Number of coefficients does not match with the design matrix.");

    apply(covX, 1, function(x) {sum(x * regCoeff)});
}


get.covmat <- function(StDevCovar, corrCovar) {
    n.x      <- length(StDevCovar);
    Vars     <- StDevCovar*StDevCovar;
    CovarMat <- matrix(NA, n.x, n.x);
    for (i in 1:n.x) {
        CovarMat[i,i] <- Vars[i];
        for (j in i:n.x) {
            if (j == i) {
                CovarMat[i,i] <- Vars[i];
                next;
            }
            CovarMat[i, j] <- corrCovar*StDevCovar[i]*StDevCovar[j];
            CovarMat[j, i] <- CovarMat[i, j];
        }
    }

    CovarMat
}

## cut covariates into categories
get.cov.cat <- function(covX, breaks = NULL) {
    f.cut <- function(x, bs) {
        if (is.null(bs))
            return(x);

        bs  <- sort(unique(c(-Inf, bs, Inf)));
        rst <- as.numeric(cut(x, breaks = bs)) - 1;
        factor(rst);
    }

    if (is.null(breaks))
        return(covX);

    if (is.numeric(breaks)) {
        rst <- apply(covX, 1, f.cut, breaks);
        rst <- t(rst);
    } else if (is.list(breaks)) {
        rst <- covX;
        for (i in 1:min(ncol(covX), length(breaks))) {
            rst[,i] <- f.cut(covX[,i], breaks[[i]]);
        }
    }

    rst
}
