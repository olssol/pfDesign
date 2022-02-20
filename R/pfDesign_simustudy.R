#' Simulate patients for conduction simulation studies
#'
#' @param simu_setting A list of simulation parameters for the following
#'     functions \code{pdSimuPts}, \code{pdGetInt}
#' @param h1 Simulate patients under null (h1 = FALSE) or alternative hypothesis
#'     (h1 = TRUE)
#'
#'
#' @examples
#'     simu_setting <- list(n_pat      = 100,
#'                          mu_0       = -1.5,
#'                          mu_0_trt   = -0.4,
#'                          fmla_t     = as.formula("~ -1 + time"),
#'                          coeff_t    = 0,
#'                          muCov      = rep(0,   5)
#'                          sdCov      = rep(0.5, 5)
#'                          corCov     = 0.1,
#'                          coeff_x    = rep(0.5, 5)
#'                          t_end      = 3,
#'                          t_start    = 0,
#'                          t_hist     = 2,
#'                          intervals  = 2,
#'                          seed       = 1000)
#'
#'    simu_patients(simu_setting)
#'
#' @export
#'
pd_ss_patients <- function(simu_setting, h1 = FALSE) {

    sum_truth <- function(dta) {
        mi <- max(dta$interval)
        dta %>%
            group_by(interval, itime, ilabel) %>%
            summarise(n  = n(),
                      mt = mean(time),
                      mx = mean(xbeta),
                      mg = mean(tgamma),
                      sy = sum(y),
                      y  = mean(y))
    }

    ## under alternative for treatment arm
    if (h1) {
        print("Under H1")
        simu_setting$mu_0 <- simu_setting$mu_0_trt
    }

    all_Pts   <- do.call(pdSimuPts, simu_setting)
    intervals <- simu_setting$intervals

    rst <- list()
    for (i in seq_len(length(intervals))) {
        Intv <- do.call(pdGetInt,
                        c(simu_setting,
                          list(vec_t      = all_Pts$time,
                               n_interval = intervals[i])
                          ))

        Pts  <- cbind(all_Pts, Intv)
        Sum  <- sum_truth(Pts)

        rst[[i]] <- list(pts        = Pts,
                         summary    = Sum,
                         n_interval = intervals[i])
    }

    rst
}

#'  Simulation Truth
#'
#' @inherit pd_ss_patients
#'
#' @param ntruth number of patients to simulate for summarizing the truth
#'
#'
#' @export
#'
pd_ss_truth <- function(simu_setting, ntruth =1000000, ...) {
    true_setting       <- simu_setting
    true_setting$n_pat <- ntruth
    rst_truth          <- pd_ss_patients(true_setting, ...)

    rst_truth
}

#' Draw inference for Epsilon = 0 and 1
#'
#' @param dta_sum summary of the simulated patients, a return from
#'     \code{pd_ss_patients}
#'
#'
#' @export
#'
pd_ss_get_reference <- function(dta_sum, labels = c("Eps = 0", "Eps = 1")) {
    ## Typical Posterior: Interval Only
    rst_theta_int <- NULL
    for (i in seq_len(nrow(dta_sum))) {
        cur_n      <- as.numeric(dta_sum[i, "n"])
        cur_sy     <- as.numeric(dta_sum[i, "sy"])
        cur_theta  <- rbeta(4000,
                            cur_sy + 0.25,
                            cur_n - cur_sy + 0.25)
        rst_theta_int <- rbind(rst_theta_int,
                               data.frame(theta    = cur_theta,
                                          interval = i))
    }
    rst_theta_int$grp     <- labels[1]
    rst_theta_int$epsilon <- 0

    ## Typical Posterior: Cumulative Data
    rst_theta_cum <- NULL
    for (i in seq_len(nrow(dta_sum))) {
        cur_n_cum  <- as.numeric(sum(dta_sum[1:i, "n"]))
        cur_sy_cum <- as.numeric(sum(dta_sum[1:i, "sy"]))
        cur_theta_cum <- rbeta(4000,
                               cur_sy_cum + 0.25,
                               cur_n_cum - cur_sy_cum + 0.25)
        rst_theta_cum <- rbind(rst_theta_cum,
                               data.frame(theta    = cur_theta_cum,
                                          interval = i))
    }
    rst_theta_cum$grp     <- labels[2]
    rst_theta_cum$epsilon <- 1

    ## combine and return
    rst_theta         <- rbind(rst_theta_int, rst_theta_cum)
    rst_theta$type    <- "Reference"

    rownames(rst_theta) <- NULL
    rst_theta
}

#' Draw inference by filtering
#'
#'  @param dta_pts Simulated patients, a return from \code{{pd_ss_patients}
#'  @param v_eps Vector of epsilon values
#'
#' @export
#'
pd_ss_get_filtering <- function(dta_pts, v_eps) {
    rst_filter <- NULL;
    for (eps in v_eps) {
        cur_theta      <- pdFilter(dta_pts$y, dta_pts$interval,
                                   epsilons = eps,
                                   nsmps = 100000,
                                   mix_ab = c(1, 1))
        cur_theta$grp <- paste("Eps = ", eps, sep = "");
        rst_filter    <- rbind(rst_filter, cur_theta)
    }

    rownames(rst_filter) <- NULL
    rst_filter
}

#' Draw posterior distribution through Beta
#'
#'
#' @export
#'
pd_ss_get_betapost <- function(trt_pts, npost = 3000) {
    pts <- trt_pts[[1]]$pts %>%
        dplyr::filter(ilabel == "Concurrent")

    y   <- pts$y
    rst <- rbeta(npost, sum(y) + 0.5, sum(1 - y) + 0.5)
    rst
}



#' Conduct logistic regression with a parametric model
#'
#' @inherit pd_ss_get_filtering
#'
#' @param nx number of covariates to put in the model
#'
#' @export
#'
pd_ss_get_para <- function(dta_pts, nx = 0) {
    fml <- "y ~ time"
    if (nx > 0) {
        fml <- paste(c(fml, paste("V", 1:nx, sep = "")), collapse = "+")
    }

    rst_par <- pdLogist(dta = dta_pts, fml = formula(fml))
    rst_par
}

#' Summarize simulation study results
#'
#'
#'
#' @export
#'
pd_ss_summarize <- function(rst_truth, rst_theta, rst_filter, rst_par,
                            rst_par0, rst_trt, n_current,
                            cur_set = NULL) {

    ## treatment vs. control p-value
    get_pval <- function(est1, var1, est2, var2) {
        z    <- est2 - est1
        z    <- z / sqrt(var1 + var2)
        pz   <- pnorm(abs(z))
        pval <- 2 * min(pz, 1 - pz)
        p_g0 <- 1 - pnorm(0, z)

        c(pval, p_g0)
    }

    ## treatment vs. control eff > 0 posterior probability
    get_eff <- function(post1, post2) {
        n        <- max(length(post1), length(post2))
        post1    <- sample(post1, size = n, replace = T)
        post2    <- sample(post2, size = n, replace = T)
        post_eff <- post2 - post1
        p_g0     <- mean(post_eff > 0)
        pval     <- get_pval(mean(post1), var(post1),
                             mean(post2), var(post2))

        c(pval, p_g0)
    }

    f_theta <- function(rst, g) {
        int <- max(rst$interval)
        rst <- rst %>%
            dplyr::filter(grp      == g   &
                          interval == int &
                          type     == "Posterior")

        rst$theta
    }

    f_filter_ess <- function(rst, g, theta_eps0) {
        int   <- max(rst$interval)
        theta <- f_theta(rst, g)
        ess   <- pdESS(theta, theta_eps0, n_current)

        data.frame(grp        = g,
                   type       = "Posterior",
                   interval   = int,
                   ess        = ess[1],
                   ess_tot    = ess[2])
    }

    f_filter_pow <- function(rst, g) {
        int   <- max(rst$interval)
        theta <- f_theta(rst, g)

        ## type I
        type1 <- get_eff(theta, rst_trt[[1]])

        ## power
        power <- get_eff(theta, rst_trt[[2]])

        data.frame(grp        = g,
                   type       = "Posterior",
                   interval   = int,
                   type1_pval = type1[1],
                   type1_prob = type1[2],
                   power_pval = power[1],
                   power_prob = power[2])

    }

    f_par_pow <- function(rst, grp) {
        interval <- rst[nrow(rst), "interval"]
        est1     <- rst[nrow(rst), "mean"]
        var1     <- rst[nrow(rst), "variance"]

        ## type I
        est2  <- mean(rst_trt[[1]])
        var2  <- var(rst_trt[[1]])
        type1 <- get_pval(est1, var1, est2, var2)

        browser()

        ## power
        est2  <- mean(rst_trt[[2]])
        var2  <- var(rst_trt[[2]])
        power <- get_pval(est1, var1, est2, var2)

        data.frame(grp  = grp,
                   type = grp,
                   interval   = interval,
                   type1_pval = type1[1],
                   type1_prob = type1[2],
                   power_pval = power[1],
                   power_prob = power[2])
    }

    f_par <- function(rst, grp) {
        data.frame(grp      = grp,
                   type     = grp,
                   interval = rst[, "interval"],
                   mean     = rst[, "mean"],
                   variance = rst[, "variance"],
                   lb       = 0,
                   ub       = 1)
    }

    ## ESS
    rst_ess    <- NULL
    theta_eps0 <- f_theta(rst_filter, "Eps = 0")
    for (g in unique(rst_filter$grp)) {
        rst_ess <- rbind(rst_ess,
                         f_filter_ess(rst_filter,
                                      g,
                                      theta_eps0))
    }


    ## power and type I error
    rst_power <- f_par_pow(rst_par,  "ParaWX") %>%
        rbind(f_par_pow(rst_par0, "ParaWoX"))

    for (g in unique(rst_filter$grp)) {
        rst_power <- rst_power %>%
            rbind(f_filter_pow(rst_filter, g))
    }

    ## control arm estimation
    rst_bayesian <- rbind(rst_theta, rst_filter)
    cur_summary  <- rst_bayesian %>%
        group_by(grp, type, interval) %>%
        summarise(mean     = mean(theta),
                  variance = var(theta),
                  lb       = quantile(theta, 0.025),
                  ub       = quantile(theta, 0.975)) %>%
        ## dplyr::filter(!(type %in% c("Reference"))) %>%
        ungroup(grp, type) %>%
        rbind(f_par(rst_par,  "ParaWX")) %>%
        rbind(f_par(rst_par0, "ParaWoX")) %>%
        mutate(grp  = factor(grp),
               type = factor(type)) %>%
        left_join(rst_truth  %>%
                  mutate(truth = y) %>%
                  select(interval, truth),
                  by = "interval") %>%
        ## power and type I
        left_join(rst_power) %>%
        ## ess
        left_join(rst_ess)

    ## return
    rst <- cur_summary
    if (!is.null(cur_set)) {
        rst <- cbind(rst,
                     cur_set[rep(1,
                                 nrow(cur_summary)),
                             ])
    }

    rownames(rst) <- NULL
    rst
}

#' Get coefficients for a polynomial model to mimic a sigmoid curve
#'
#'
#' @export
#'
ps_ss_get_par <- function(t_start, t_end, fml, L = 0.4, k = 1, t0 = 2) {
    fml    <- paste(c("y ", as.character(fml)), collapse = "")
    time   <- seq(t_start, t_end, length = 1000)
    true_y <- L/(1+exp(-k * (time - t0)))

    e_lm   <- lm(as.formula(fml), data.frame(time = time, y = true_y))
    pred_y <- predict(e_lm)

    plot(time, true_y, col = "green")
    lines(time, pred_y, col = 'red')

    print(e_lm$coefficients)
}


#' Combine simulation study results
#'
#'
#'
#' @export
#'
pd_ss_combine <- function(prefix, cmb_reps = 1:50, f_rst = "rst.Rdata") {
    if (file.exists(f_rst)) {
        print(paste(f_rst, " exists...", sep = ""));
        return(NULL);
    }

    ##combine
    all_rst <- NULL;
    for (cr in seq_len(length(cmb_reps))) {
        r     <- cmb_reps[cr]
        f_fit <- paste(prefix, "_", r, ".Rdata", sep = "");
        if (!file.exists(f_fit)) {
            print(paste(f_fit, " does not exist...", sep = ""))
            next
        } else {
            load(f_fit)
        }

        all_rst <- rbind(all_rst, simu_rst)
    }

    save(all_rst, file = f_rst)
}

#' Generate plot for posterior distributions
#'
#'
#' @export
#'
pd_ss_plt_bayes <- function(rst_bayesian, rst_summary, n_interval = 4) {
    ## rst_truth <- rst_summary %>%
    ##   select(interval, truth) %>%
    ##   dplyr::filter(interval < 5) %>%
    ##   mutate(interval = paste("Interval~", interval)) %>%
    ##   distinct()

    ss <- rst_bayesian$type
    ss[which(rst_bayesian$type == "Mixed")] <- "TEA"
    ss[which(rst_bayesian$type == "Reference"
             & rst_bayesian$grp == "Eps = 0")] <- "Reference (Ignoring)"
    ss[which(rst_bayesian$type == "Reference" &
             rst_bayesian$grp == "Eps = 1")] <- "Reference (Pooling)"
    rst_bayesian$type  <- ss

    rst_plot <- rst_bayesian %>%
        dplyr::filter(interval <= n_interval) %>%
        mutate(interval = paste("Interval~", interval),
               grp      = factor(grp,
                                 labels = paste("epsilon ==",
                                                c(0, 0.25, 0.50, 0.75, 1),
                                                sep = "")))

    plt_1 <-  ggplot(data = rst_plot, aes(x = theta)) +
        stat_density(aes(group=type, color = type),
                     position="identity", geom = "line", adjust = 2) +
        ## geom_vline(data = rst_truth, aes(xintercept = truth)) +
        theme_bw() +
        theme(strip.background = element_blank(),
              legend.title= element_blank(),
              legend.position = "bottom") +
        scale_x_continuous(limits = c(0.05, 0.95)) +
        labs(x = expression(theta), y = "Density") +
        facet_grid(grp ~ interval, labeller = label_parsed)

    plt_1
}

#' Plot observed response rates
#'
#'
#'
#' @export
pd_ss_plt_truth <- function(dtaSum, fname = "plt_interval.pdf",
                            ylim = c(0.25, 0.36), xlim = c(0, 3)) {

    dtaSum <- dtaSum %>%
      mutate(Scenario = paste("Scenario", as.roman(Scenario)),
             NInt     = paste("#Intervals = ", NInt))

    plt  <-  ggplot(data = dtaSum, aes(x = itime, y = y)) +
      geom_line(lty = 2, col = "gray60") +
      geom_point(aes(shape = ilabel)) +
      theme_bw() +
      ##theme(axis.text.x = element_text(color = "black", size=9, angle=0,
      ##                                 vjust=0, hjust=0.7)) +
      theme(legend.position = "bottom",
            panel.grid = element_blank(),
            legend.title = element_blank(),
            strip.background = element_blank()) +
      scale_y_continuous(limits = ylim) +
      scale_x_continuous(limits = xlim) +
      labs(y = "Average Response Rate", x = "Time") +
      facet_grid(NInt ~ Scenario)

    d_text <- rbind(data.frame(x = 1, y = 0.1, label = "Non-concurrent"),
                    data.frame(x = 2.5, y = 0.1, label = "Concurrent"))

    plt <- plt +
        geom_vline(xintercept = 2, lty = 3, color = "red") +
        geom_text(
            data = d_text,
            aes(x = x, y = y, label = label),
            color = "gray30",
            angle = 0,
            size = 2
        )

    plt
}


#' Entire simulation process
#'
#'
#'
#' @export
#'
pd_ss_simu_all <- function(simu_setting, rst_truth, seed = NULL, ...) {

    if (!is.null(seed)) {
        old_seed <- set.seed(seed)
    }

    ##-----------------------------------------------------------------
    ##             SIMU PATIENTS
    ##-----------------------------------------------------------------
    ## control
    lst_dta_pts            <- pd_ss_patients(simu_setting)
    ## treatment arm
    simu_setting_trt       <- simu_setting
    simu_setting_trt$n_pat <- g_n_pat_trt
    ## treatment under H0
    lst_dta_pts_trt_h0     <- pd_ss_patients(simu_setting_trt)
    ## treatment under H1
    simu_setting_trt$mu_0  <- simu_setting$mu_0_trt
    lst_dta_pts_trt_h1     <- pd_ss_patients(simu_setting_trt)

    ##-----------------------------------------------------------------
    ##             ESTIMATION: TREATMENT
    ##-----------------------------------------------------------------
    rst_trt <- list(pd_ss_get_betapost(lst_dta_pts_trt_h0),
                    pd_ss_get_betapost(lst_dta_pts_trt_h1))

    ##-----------------------------------------------------------------
    ##             ESTIMATION: Control and Comparison
    ##-----------------------------------------------------------------
    simu_rst <- NULL
    for (i in seq_len(length(lst_dta_pts))) {
        dta_sum     <- lst_dta_pts[[i]]$summary
        dta_pts     <- lst_dta_pts[[i]]$pts
        n_current   <- length(which(dta_pts$ilabel == "Concurrent"))
        n_interval  <- lst_dta_pts[[i]]$n_interval

        ## reference
        rst_theta  <- pd_ss_get_reference(dta_sum)

        ## filtering
        rst_filter <- pd_ss_get_filtering(dta_pts, g_eps)

        ## parametric with x
        rst_par    <- pd_ss_get_para(dta_pts, g_nx)

        ## parametric without x
        rst_par0   <- pd_ss_get_para(dta_pts, 0)

        ##-----------------------------------------------------------------
        ##             SUMMARY RESULTS
        ##-----------------------------------------------------------------
        rst_summary <- pd_ss_summarize(rst_truth[[i]]$summary,
                                       rst_theta, rst_filter,
                                       rst_par, rst_par0, rst_trt,
                                       n_current, ...)

        rst_summary <- rst_summary %>%
            mutate(inx_pt = i) %>%
            dplyr::filter(interval == n_interval + 1)

        ##-----------------------------------------------------------------
        ##             APPEND
        ##-----------------------------------------------------------------
        simu_rst <- rbind(simu_rst, rst_summary)
    }

    ##-----------------------------------------------------------------
    ##             RETURN
    ##-----------------------------------------------------------------
    if (!is.null(seed)) {
        set.seed(old_seed)
    }

    simu_rst
}
