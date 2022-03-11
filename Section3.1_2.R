# ------------------------------------------------------------------------------------------------------------------ 
# Section 3. Simulation Study - 3.1. Performance of Simultaneous Confidence Bands - second example
# ------------------------------------------------------------------------------------------------------------------

# ------------ load the packages -------------- 
rm(list=ls())
#dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)  # create personal library
#.libPaths(Sys.getenv("R_LIBS_USER"))  # add to the path
if (!require(parallel)) {
  install.packages("parallel")
  require(parallel) # this package is needed for the \code{makeCluster} and \code{stopCluster} functions below for cluster computing
}
if (!require(foreach)) {
  install.packages("foreach")
  require(foreach) # we use version 1.5.1 of this package, which is needed for the \code{foreach} function below for cluster computing
}
if (!require(doParallel)) {
  install.packages("doParallel")
  require(doParallel) # we use version 1.0.16 of this package, which is needed for the \code{registerDoParallel} function below for cluster computing
}
if (!require(devtools)) {
  install.packages("devtools") # we use version 2.4.2 of this package, which is needed for installing the fdEL and fregion packages below
  require(devtools)
}
if (!require(fdEL)) {
  install_github("news11/fdEL", force = TRUE)
  require(fdEL) # we use version 0.0.0.9000 of this package, which is needed for implementing the EL procedures, the Cao and Cao2 confidence bands, and the Geo and Fmax tests
} 
if (!require(SCBmeanfd)) {
  install.packages("SCBmeanfd")
  require(SCBmeanfd) # we use version 1.2.2 of this package, which is needed for implementing the MFD confidence band
}
if (!require(KernSmooth)) {
  install.packages("KernSmooth")
  require(KernSmooth) # we use version 2.23-18 of this package, which is needed for the locpoly function inside the function scb.mean00
}
if (!require(fregion)) {
  install_github("hpchoi/fregion", force = TRUE)
  require(fregion) # we use version 0.0934 of this package, which is needed for implementing the Geo confidence band
}

# ------------ user input -------------- #
dir_path = "C:\\R_HW\\paper_fnl_codes_20220111\\source" # where the source files 'scb.mean00.R' and 'Geoband.R' are; please manually download them from "https://github.com/news11/paper_fdEL" 
dir_path2=paste("C:\\R_HW\\fnl_1CB_20180602file_outs_Xtinteger_20210703",sep="") # where the resulting files can be saved
parameters = list()
# Results in Table 1, fifth column correspond to setting n_subject = 100, by_agrid = 4, sigma = 1.5 below
# Results in Table 1, sixth column correspond to setting n_subject = 200, by_agrid = 2, sigma = 1.5 below
# Results in Table 1, seventh column correspond to setting n_subject = 100, by_agrid = 4, sigma = 2 below
# Results in Table 1, eighth column correspond to setting n_subject = 200, by_agrid = 2, sigma = 2 below
parameters$no_cores <- 50 # the number of cores for parallel computing, preferably < parallel::detectCores(); changing it can still reproduce the results in case 2 of table 1, as long as parameters$no_cores * parameters$nrep_sub = 1000
parameters$nrep_sub = 20 # number of datasets to be generated per parallel task; changing it can still reproduce the results in case 2 of table 1, as long as parameters$no_cores * parameters$nrep_sub = 1000
parameters$n_subject = 100 # number of subjects
parameters$nboot = 1000 # number of bootstrap samples 
parameters$alpha_vec = 0.05 # a number representing the significance level of interest
parameters$tau = 1 # [0, tau): the time domain of the \eqn{X(t)} process
parameters$sigma = 1.5 # a number with values 1.5 or 2, representing \eqn{\nu_L} in simulating the log-normally distributed \eqn{\epsilon}
parameters$by_agrid = 4  # a number with values 2 (when \code{parameters$n_subject == 200}) or 4 (when \code{parameters$n_subject == 100}), representing the mesh of \eqn{{\bf G}_n} (defined in Section 2.2, from the smallest to the largest)
# 20220304 parameters$tol_resultEEs = 10 ^ (-4), parameters$n_lamb_grid = 1000, parameters$m1m2lbub_digit = 6, parameters$var_option = 0
parameters$right_endpt = c(0.05) # the value(s) of \eqn{z} in Supplement Section 5.1.
parameters$by_agrid_eval = 1 # the mesh of the set of activity levels at which we want to evaluate the confidence band
parameters$c1 = (100000) ^ 2 # the constant 316 in simulating \eqn{X(t)}, raised to the fourth power 
parameters$c2 = (1 / 1000) ^ 2 # the constant added in simulating \eqn{epsilon}
parameters$gridsol = parameters$tau / (10 * parameters$n_subject) # the mesh of \eqn{{\bf H}_n} (defined in Section 3.1)

# ------------ define functions to be used -------------- #
setwd(dir_path) # for sourcing files in the following lines
source('scb.mean00.R', local = TRUE)
source('Geoband.R', local = TRUE)
floor_n = function(x, n = 1) {
  # Extending the floor function so that for the given input number x and positive integer n,
  # the function returns the number b / n for some interger b, such that b / n <= x < (b + 1) / n
  # The round function is to eliminate a small numerical difference from the actual value of x * n.
  # Args:
  #   x : a real number
  #   n: a positive integer n
  # Returns:
  #   the number b / n for some interger b, such that b / n <= x < (b + 1) / n
  floor(round(x * n, 2)) / n
}
coveragefn = function(parameters, n_subject, sigma, split) {
  # Simulate data from the second example in Section 3.1 and computes the empirical coverage of EL, EP and NS confidence bands
  # Args:
  #   parameters: a list of parameters from the ``user input'' above
  #   n_subject: the number of subjects 
  #   split: number of parallel tasks
  # Returns: 
  #   a list containing the following:
  #   CB_out_all: a list whose \code{nrepi}-th element contains some of the outputs from the functions for computing the EL, EP, NS, Cao, Cao2, MFDbs, Geo, naivet and MFD confidence bands based on the \code{nrepi}-th simulated dataset
  #   coverages: the empirical coverage of EL, EP, NS, Cao, Cao2, MFDbs, Geo, naivet and MFD confidence bands
  t_ncol = parameters$tau / parameters$gridsol  
  right_endpt_len = length(parameters$right_endpt)
  # define output
  EL_CB_Nairrej_all = matrix(0, nrow = parameters$nrep_sub, ncol = right_endpt_len)
  EP_CB_Nairrej_all = matrix(0, nrow = parameters$nrep_sub, ncol = right_endpt_len)
  HW_CBrej_all = matrix(0, nrow = parameters$nrep_sub, ncol = right_endpt_len)
  CB_out_all = list()
  Cao_CB0rej_all = matrix(0, nrow = parameters$nrep_sub, ncol = right_endpt_len)
  Cao_CBrej_all = matrix(0, nrow = parameters$nrep_sub, ncol = right_endpt_len)
  Cao_CB2rej_all = matrix(0, nrow = parameters$nrep_sub, ncol = right_endpt_len)
  Geo_CBrej_all = array(0, c(parameters$nrep_sub, 3, right_endpt_len)) 
  MFD_CBrej_all = matrix(0, nrow = parameters$nrep_sub, ncol = right_endpt_len)
  # for loop generating \code{parameters$nrep_sub} datasets and computing the confidence bands accordingly
  for (nrepi in 1:parameters$nrep_sub) {
    set.seed(parameters$nrep_sub * (split - 1) + nrepi, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection") 
    # simulate \eqn{X(t)} in the second example in Section 3.1
    t = t(sapply(1:n_subject, FUN = function (i) {
      sample(seq(0, parameters$tau - parameters$gridsol, by = parameters$gridsol), t_ncol, replace = FALSE)
    }))
    eps = rlnorm(n_subject, meanlog = 0, sdlog = sigma) + parameters$c2   
    eps_tmat = matrix(eps, nrow = n_subject, ncol = t_ncol)
    Unif_rv = runif(n_subject, 0, parameters$tau)
    Unif_rv_tmat = matrix(Unif_rv, nrow = n_subject, ncol = t_ncol)
    Ind_t_g_U = (t > Unif_rv_tmat)
    Xt = floor_n((parameters$c1 * t / eps_tmat) ^ 0.25) * Ind_t_g_U
    as_all = seq(0, (parameters$tau * parameters$c1 / parameters$c2) ^ 0.25, by = 1) 
    as = seq(range(as_all)[1], range(as_all)[2], by = parameters$by_agrid)
    as_eval = seq(range(as_all)[1], range(as_all)[2], by = parameters$by_agrid_eval)
    n_agrid_eval = length(as_eval)
    c_a_vec = t(sapply(1:n_agrid_eval, FUN = function (i) {
      if (floor_n(as_eval[i]) == as_eval[i] & as_eval[i] < range(as_all)[2]) {  
        return(as_eval[i] + 1)  
      } else {
        return(ceiling(as_eval[i]))
      }
    }))
    # compute the true mean occupation time \eqn{E\{L(a)\} in the second example in Section 3.1; the computation follows the procedure described in the second paragraph of Supplement Section 12
    mu_1 = 0 + log((c_a_vec ^ (1 / 0.25)) / parameters$c1)
    ba_1p_c2_c1 = ((c_a_vec ^ (1 / 0.25)) * parameters$c2 / parameters$c1)
    H_n = c(seq(0, parameters$tau - parameters$gridsol, by = parameters$gridsol), Inf)
    P_ba2_ge = t(sapply(1:length(H_n), FUN = function (i) {  # \eqn{P(c_a^H \ge t), t \in {\bf H}_n} in the second paragraph of Supplement Section 12
      if (i == 1) {
        return(rep(1, times = n_agrid_eval))
      } else {
        xm1_in_H_n = H_n[i - 1]
        x_1 = xm1_in_H_n - ba_1p_c2_c1 # length n_agrid_eval vec
        x1_mean_sd = (log(x_1 * (x_1 >= 0)) - mu_1) / sigma 
        return(1 - pnorm(x1_mean_sd))
      }
    }))
    P_ba3_ge = t(sapply(1:length(H_n), FUN = function (i) { # \eqn{P(U_H \ge t), t \in {\bf H}_n} in the second paragraph of Supplement Section 12
      if (i == 1) {
        return(rep(1, times = n_agrid_eval))
      } else {
        xm1_in_H_n = H_n[i - 1]
        return(rep(1 - xm1_in_H_n / parameters$tau, times = n_agrid_eval))
      }
    }))
    P_ba2_g = rbind(P_ba2_ge[-1, ], 0) # \eqn{P(c_a^H > t), t \in {\bf H}_n} in the second paragraph of Supplement Section 12
    P_ba3_g = rbind(P_ba3_ge[-1, ], 0) # \eqn{P(U_H > t), t \in {\bf H}_n} in the second paragraph of Supplement Section 12
    P_ba2_l = 1 - P_ba2_ge # \eqn{P(c_a^H < t), t \in {\bf H}_n} in the second paragraph of Supplement Section 12
    P_ba3_l = 1 - P_ba3_ge # \eqn{P(U_H < t), t \in {\bf H}_n} in the second paragraph of Supplement Section 12
    P_ba2_le = 1 - P_ba2_g # \eqn{P(c_a^H \le t), t \in {\bf H}_n} in the second paragraph of Supplement Section 12
    P_ba3_le = 1 - P_ba3_g # \eqn{P(U_H \le t), t \in {\bf H}_n} in the second paragraph of Supplement Section 12
    P_A_le = P_ba2_le * P_ba3_le # \eqn{P(A_H \le t), t \in {\bf H}_n} in the second paragraph of Supplement Section 12
    P_A_l = P_ba2_l * P_ba3_l # \eqn{P(A_H < t), t \in {\bf H}_n} in the second paragraph of Supplement Section 12
    P_A = P_A_le - P_A_l # \eqn{P(A_H = t), t \in {\bf H}_n} in the second paragraph of Supplement Section 12
    rm(mu_1, ba_1p_c2_c1, P_ba2_ge, P_ba3_ge, P_ba2_g, P_ba3_g, P_ba2_le, P_ba3_le, P_ba2_l, P_ba3_l, P_A_le, P_A_l)
    mu_a = apply(matrix(pmax(parameters$tau - H_n, 0), byrow = FALSE, nrow = length(H_n), ncol = n_agrid_eval) * P_A, 2, sum) # \eqn{E\{L(a)\} in the second paragraph of Supplement Section 12
    rm(H_n, P_A)
    # compute \eqn{L(a)} in the second example in Section 3.1
    La = t(sapply(1:n_subject, FUN = function (i) {  
      fdEL::mu_a_hat(Xt[i, ], grid_widths = rep(parameters$gridsol, t_ncol), as) 
    }))
    mu_a0 = apply(La, 2, mean)
    b_a_vec = t(sapply(1:n_agrid_eval, FUN = function (i) {
      min(which(as >= as_eval[i]))
    }))
    right_ind_keep = sapply(parameters$right_endpt, function(i) {
      if (sum(mu_a0[b_a_vec] / parameters$tau > i) != 0) {  
        max(which(mu_a0[b_a_vec] / parameters$tau > i))
      } else {1}
    })  
    right_ind_keep_as = sapply(parameters$right_endpt, function(i) {
      if (sum(mu_a0 / parameters$tau > i) != 0) {
        max(which(mu_a0 / parameters$tau > i))
      } else {1}
    })  
    # computing the EL confidence band
    CB_out = fdEL::elfband(Ta = La, as = as, as_eval = as_eval, n_boot = parameters$nboot, mu_a = mu_a, alpha = parameters$alpha_vec, as_right_ind_keep0 = right_ind_keep, as_right_ind_keep_as0 = right_ind_keep_as)
    EL_CB_Nairrej_all[nrepi, ] = CB_out$EL_CB_Nairrej
    EP_CB_Nairrej_all[nrepi, ] = CB_out$EP_CB_Nairrej
    HW_CBrej_all[nrepi, ] = CB_out$HW_CBrej
    # computing the MFDbs confidence band mentioned in Section 3
    GeoCB_Bs <- fmean(Y = t(La), as = as, as_eval = as_eval, alpha = parameters$alpha_vec[1], btypes = c("Bs"), as_right_ind_keep_as0 = right_ind_keep_as, right_endpt = parameters$right_endpt)
    for (as_right_ind in 1:right_endpt_len) {
      # not used: computing the Cao confidence band before projecting the initially smoothed covariance estimates onto the space of non-negative definite matrices    
      an.error.occured <- FALSE
      a.warning.occured <- FALSE
      possibleError_Cao0 <- tryCatch(
        CaoCB0 <- fdEL::Caoband(Y = La[, 1:right_ind_keep_as[as_right_ind]], as = as[1:right_ind_keep_as[as_right_ind]], checkpts = as_eval[1:right_ind_keep[as_right_ind]], b_M = parameters$nboot, mean_true = mu_a[1:right_ind_keep[as_right_ind]], alpha = parameters$alpha_vec[1], cov_nnd = 0),
        error = function(e) {an.error.occured <<- TRUE},
        warning = function(w) {a.warning.occured <<- TRUE}
      )  
      if (sum(an.error.occured + a.warning.occured) > 0) {
        CB_out$CaoCB0[[as_right_ind]] = matrix(-1, ncol = 2, nrow = right_ind_keep[as_right_ind]) 
      } else {
        CB_out$CaoCB0[[as_right_ind]] = CaoCB0$confidence_bounds
      } 
      CB_out$CaoCB0_Nair[[as_right_ind]] = rbind(CB_out$CaoCB0[[as_right_ind]], matrix(-1, ncol = 2, nrow = n_agrid_eval - right_ind_keep[as_right_ind]))
      CB_out$CaoCB0_Nair[[as_right_ind]][-(1:right_ind_keep[as_right_ind]), 1] = (CB_out$CaoCB0[[as_right_ind]])[right_ind_keep[as_right_ind], 1]
      CaoCB0_as_right = (CB_out$CaoCB0_Nair[[as_right_ind]])[right_ind_keep[as_right_ind], 2]
      CB_out$CaoCB0_Nair[[as_right_ind]][-(1:right_ind_keep[as_right_ind]), 2] = CaoCB0_as_right * (CaoCB0_as_right < 0)
      Cao_CB0rej_all[nrepi, as_right_ind] = (sum(((CB_out$CaoCB0_Nair[[as_right_ind]][, 1] >= mu_a) * (CB_out$CaoCB0_Nair[[as_right_ind]][, 2] <= mu_a)) == 0) >= 1)  
      # computing the Cao confidence band mentioned in Section 3    
      CaoCB <- fdEL::Caoband(Y = La[, 1:right_ind_keep_as[as_right_ind]], as = as[1:right_ind_keep_as[as_right_ind]], checkpts = as_eval[1:right_ind_keep[as_right_ind]], b_M = parameters$nboot, mean_true = mu_a[1:right_ind_keep[as_right_ind]], alpha = parameters$alpha_vec[1])
      CB_out$CaoCB[[as_right_ind]] = CaoCB$confidence_bounds
      CB_out$CaoCB_Nair[[as_right_ind]] = rbind(CB_out$CaoCB[[as_right_ind]] ,matrix(-1, ncol = 2, nrow = n_agrid_eval - right_ind_keep[as_right_ind]))
      CB_out$CaoCB_Nair[[as_right_ind]][-(1:right_ind_keep[as_right_ind]), 1] = (CB_out$CaoCB[[as_right_ind]])[right_ind_keep[as_right_ind], 1]
      CaoCB_as_right = (CB_out$CaoCB_Nair[[as_right_ind]])[right_ind_keep[as_right_ind], 2]
      CB_out$CaoCB_Nair[[as_right_ind]][-(1:right_ind_keep[as_right_ind]), 2] = CaoCB_as_right * (CaoCB_as_right < 0)
      Cao_CBrej_all[nrepi, as_right_ind] = (sum(((CB_out$CaoCB_Nair[[as_right_ind]][, 1] >= mu_a) * (CB_out$CaoCB_Nair[[as_right_ind]][, 2] <= mu_a)) == 0) >= 1)  
      # computing the Cao2 confidence band mentioned in Section 3    
      CaoCB2 <- fdEL::Caoband(Y = La[, 1:right_ind_keep_as[as_right_ind]], as = as[1:right_ind_keep_as[as_right_ind]], checkpts = as_eval[1:right_ind_keep[as_right_ind]], b_M = parameters$nboot, mean_true = mu_a[1:right_ind_keep[as_right_ind]], alpha = parameters$alpha_vec[1], N_m = right_ind_keep_as[as_right_ind] - 2) 
      CB_out$CaoCB2[[as_right_ind]]= CaoCB2$confidence_bounds
      CB_out$CaoCB2_Nair[[as_right_ind]] = rbind(CB_out$CaoCB2[[as_right_ind]], matrix(-1, ncol = 2, nrow = n_agrid_eval - right_ind_keep[as_right_ind]))
      CB_out$CaoCB2_Nair[[as_right_ind]][-(1:right_ind_keep[as_right_ind]), 1] = (CB_out$CaoCB2[[as_right_ind]])[right_ind_keep[as_right_ind], 1]
      CaoCB2_as_right = CB_out$CaoCB2_Nair[[as_right_ind]][right_ind_keep[as_right_ind], 2]
      CB_out$CaoCB2_Nair[[as_right_ind]][-(1:right_ind_keep[as_right_ind]), 2] = CaoCB2_as_right * (CaoCB2_as_right < 0)
      Cao_CB2rej_all[nrepi, as_right_ind] = (sum(((CB_out$CaoCB2_Nair[[as_right_ind]][, 1] >= mu_a) * (CB_out$CaoCB2_Nair[[as_right_ind]][, 2] <= mu_a)) == 0) >= 1)  
      # computing the MFD confidence band mentioned in Section 3
      h <- SCBmeanfd::cv.select(as[1:right_ind_keep_as[as_right_ind]], La[, 1:right_ind_keep_as[as_right_ind]], 1)  
      scbLa <- scb.mean00(as[1:right_ind_keep_as[as_right_ind]], La[, 1:right_ind_keep_as[as_right_ind]], bandwidth = h, degree = 1, gridsize = right_ind_keep[as_right_ind]) 
      CB_out$MFD_CB[[as_right_ind]] = t(scbLa$normscb)  
      CB_out$MFD_CB_Nair[[as_right_ind]] = cbind(CB_out$MFD_CB[[as_right_ind]], matrix(-1, nrow = 2, ncol = n_agrid_eval - right_ind_keep[as_right_ind]))
      CB_out$MFD_CB_Nair[[as_right_ind]][2, -(1:right_ind_keep[as_right_ind])] = (CB_out$MFD_CB[[as_right_ind]])[2, right_ind_keep[as_right_ind]]
      MFD_CB_as_right = CB_out$MFD_CB_Nair[[as_right_ind]][1, right_ind_keep[as_right_ind]]
      CB_out$MFD_CB_Nair[[as_right_ind]][1, -(1:right_ind_keep[as_right_ind])] = MFD_CB_as_right * (MFD_CB_as_right < 0)
      MFD_CBrej_all[nrepi, as_right_ind] = (sum(((CB_out$MFD_CB_Nair[[as_right_ind]][2, ] >= mu_a) * (CB_out$MFD_CB_Nair[[as_right_ind]][1, ] <= mu_a)) == 0) >= 1)  
      # computing the Geo confidence band mentioned in Section 3 (\code{btypes = "BEc"} below; \code{btypes = "naive.t"} below is not reported in Section 3 but for comparison purpose only, because it only guarantees pointwise coverage accuracy instead of simultaneous coverage accuracy)
      GeoCB_noBs <- fmean(Y = t(La[, 1:right_ind_keep_as[as_right_ind]]), as = as[1:right_ind_keep_as[as_right_ind]], as_eval = as_eval[1:right_ind_keep[as_right_ind]], alpha = parameters$alpha_vec[1], btypes = c("BEc", "naive.t"), as_right_ind_keep_as0 = right_ind_keep_as, right_endpt = parameters$right_endpt)
      CB_out$GeoCB_noBs[[as_right_ind]] = GeoCB_noBs$confidence_bounds  # columns are: BEc.u, BEc.l, naivet.u, naivet.l
      CB_out$GeoCB_Nair_noBs[[as_right_ind]] = rbind(CB_out$GeoCB_noBs[[as_right_ind]], matrix(-1, ncol = 4, nrow = n_agrid_eval - right_ind_keep[as_right_ind]))
      CB_out$GeoCB_Nair_noBs[[as_right_ind]][-(1:right_ind_keep[as_right_ind]), c(1, 3)] = matrix(CB_out$GeoCB_noBs[[as_right_ind]][right_ind_keep[as_right_ind], c(1, 3)], byrow = T, nrow = n_agrid_eval - right_ind_keep[as_right_ind], ncol = 2)
      GeoCB_as_right = CB_out$GeoCB_Nair_noBs[[as_right_ind]][right_ind_keep[as_right_ind], c(2, 4)]
      CB_out$GeoCB_Nair_noBs[[as_right_ind]][-(1:right_ind_keep[as_right_ind]), c(2, 4)] = matrix(GeoCB_as_right * (GeoCB_as_right < 0), byrow = T, nrow = n_agrid_eval - right_ind_keep[as_right_ind], ncol = 2)
      Geo_CBrej_all[nrepi, 1, as_right_ind] = (sum(((CB_out$GeoCB_Nair_noBs[[as_right_ind]][, 1] >= mu_a) * (CB_out$GeoCB_Nair_noBs[[as_right_ind]][, 2] <= mu_a)) == 0) >= 1)  
      Geo_CBrej_all[nrepi, 2, as_right_ind] = (sum(((CB_out$GeoCB_Nair_noBs[[as_right_ind]][, 3] >= mu_a) * (CB_out$GeoCB_Nair_noBs[[as_right_ind]][, 4] <= mu_a)) == 0) >= 1)  
      # continue computing the MFDbs confidence band mentioned in Section 3, based on GeoCB_Bs computed far above (before the for (as_right_ind... loop)
      CB_out$GeoCB_Bs[[as_right_ind]] = GeoCB_Bs$confidence_bounds[, ((2 * as_right_ind - 1):(2 * as_right_ind))]  
      CB_out$GeoCB_Nair_Bs[[as_right_ind]] = CB_out$GeoCB_Bs[[as_right_ind]]
      CB_out$GeoCB_Nair_Bs[[as_right_ind]][-(1:right_ind_keep[as_right_ind]), 1] = rep(CB_out$GeoCB_Bs[[as_right_ind]][right_ind_keep[as_right_ind], 1], times = n_agrid_eval - right_ind_keep[as_right_ind])
      GeoCB_as_right = CB_out$GeoCB_Nair_Bs[[as_right_ind]][right_ind_keep[as_right_ind], 2]
      CB_out$GeoCB_Nair_Bs[[as_right_ind]][-(1:right_ind_keep[as_right_ind]), 2] = rep(GeoCB_as_right * (GeoCB_as_right < 0), times = n_agrid_eval - right_ind_keep[as_right_ind])
      Geo_CBrej_all[nrepi, 3, as_right_ind] = (sum(((CB_out$GeoCB_Nair_Bs[[as_right_ind]][, 1] >= mu_a) * (CB_out$GeoCB_Nair_Bs[[as_right_ind]][, 2] <= mu_a)) == 0) >= 1)  
    }  
    CB_out_all[[nrepi]] = CB_out
  }  
  return(list(CB_out_all = CB_out_all, coverages = c(apply(1 - EL_CB_Nairrej_all, 2, mean), apply(1 - EP_CB_Nairrej_all, 2, mean), apply(1 - HW_CBrej_all, 2, mean), apply(1 - Cao_CBrej_all, 2, mean), apply(1 - Cao_CB2rej_all, 2, mean), c(t(apply(1 - Geo_CBrej_all, c(2, 3), mean))), apply(1 - MFD_CBrej_all, 2, mean))
  ))  
}  

# ------------ calculate coverage by utilizing R parallel computing -------------- #
# initiate cluster
cl = parallel::makeCluster(parameters$no_cores)
doParallel::registerDoParallel(cl)
# pass libPath to workers
parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())
# compute the results:
foreach::foreach(split = 1:parameters$no_cores,  
  .combine = c
) %dopar% {
  out = coveragefn(parameters, n_subject = parameters$n_subject, sigma = parameters$sigma, split = split)
  setwd(dir_path2)
# 20220303 still have parameters$mu = 0, parameters$t_power = 0.25, parameters$c3 = 10: save(out,  file = paste("n_", parameters$n_subject, "_split_", split, '_nrep_sub_', parameters$nrep_sub, "_mu_", parameters$mu, "_sigma_", parameters$sigma, "_gridsol_", parameters$gridsol, "_tpower_", parameters$t_power, "_by_agrid_", parameters$by_agrid, "_by_agrid_eval_", parameters$by_agrid_eval, "_c1_", parameters$c1, "_c2_", parameters$c2, "_c3_", parameters$c3, ".Rdata", sep = ""))
  save(out,  file = paste("n_", parameters$n_subject, "_split_", split, '_nrep_sub_', parameters$nrep_sub, "_sigma_", parameters$sigma, "_gridsol_", parameters$gridsol, "_by_agrid_", parameters$by_agrid, "_by_agrid_eval_", parameters$by_agrid_eval, "_c1_", parameters$c1, "_c2_", parameters$c2, ".Rdata", sep = ""))
}  
parallel::stopCluster(cl)

# ------------ collecting results from R parallel computing -------------- #
right_endpt_len = length(parameters$right_endpt)
t_ncol = parameters$tau / parameters$gridsol
Set = 1:parameters$no_cores
cov_mat = matrix(0, nrow = length(Set), ncol = 9 * right_endpt_len) 
CB_out01_mat = array(0, c(length(Set) * parameters$nrep_sub, 9, right_endpt_len))
CB_mono_mat = array(0, c(length(Set) * parameters$nrep_sub, 9, right_endpt_len))
CB_lengths_mean_mat = array(0, c(length(Set) * parameters$nrep_sub, 9, right_endpt_len))
setwd(dir_path2)
split_ind = 1
for (split in Set) {
  load(paste("n_", parameters$n_subject, "_split_", split, '_nrep_sub_', parameters$nrep_sub, "_sigma_", parameters$sigma, "_gridsol_", parameters$gridsol, "_by_agrid_", parameters$by_agrid, "_by_agrid_eval_", parameters$by_agrid_eval, "_c1_", parameters$c1, "_c2_", parameters$c2, ".Rdata", sep = ""))
  cov_mat[split_ind, ] = out$coverages
  for (sample_i in 1:parameters$nrep_sub) {
    EL_CB_Nair_out01 = 1:right_endpt_len * 0
    EP_CB_Nair_out01 = 1:right_endpt_len * 0
    HW_CB_out01 = 1:right_endpt_len * 0
    Cao_CB_out01 = 1:right_endpt_len * 0
    Cao_CB2_out01 = 1:right_endpt_len * 0
    Bs_CB_out01 = 1:right_endpt_len * 0
    BEc_CB_out01 = 1:right_endpt_len * 0
    naivet_CB_out01 = 1:right_endpt_len * 0
    MFD_CB_out01 = 1:right_endpt_len * 0
    EL_CB_Nair_length = 1:right_endpt_len * 0
    EP_CB_Nair_length = 1:right_endpt_len * 0
    HW_CB_length = 1:right_endpt_len * 0
    Cao_CB_length = 1:right_endpt_len * 0
    Cao_CB2_length = 1:right_endpt_len * 0
    Bs_CB_length = 1:right_endpt_len * 0
    BEc_CB_length = 1:right_endpt_len * 0
    naivet_CB_length = 1:right_endpt_len * 0
    MFD_CB_length = 1:right_endpt_len * 0
    EL_Nair_mono = 1:right_endpt_len * 0
    EP_Nair_mono = 1:right_endpt_len * 0
    HW_mono = 1:right_endpt_len * 0
    Cao_mono = 1:right_endpt_len * 0
    Cao2_mono = 1:right_endpt_len * 0
    Bs_mono = 1:right_endpt_len * 0
    BEc_mono = 1:right_endpt_len * 0
    naivet_mono  = 1:right_endpt_len * 0
    MFD_mono = 1:right_endpt_len * 0
    for (as_right_ind in 1:right_endpt_len) {
      EL_CB_Nair_out01[as_right_ind] = (sum((((out$CB_out_all[[sample_i]]$EL_CB_Nair[[as_right_ind]])[2, ] > parameters$tau) + ((out$CB_out_all[[sample_i]]$EL_CB_Nair[[as_right_ind]])[1, ] < 0)) > 0) >= 1)  
      EP_CB_Nair_out01[as_right_ind] = (sum((((out$CB_out_all[[sample_i]]$EP_CB_Nair[[as_right_ind]])[2, ] > parameters$tau) + ((out$CB_out_all[[sample_i]]$EP_CB_Nair[[as_right_ind]])[1, ] < 0)) > 0) >= 1)  
      HW_CB_out01[as_right_ind] = (sum((((out$CB_out_all[[sample_i]]$HW_CB_eval)[2, , as_right_ind] > parameters$tau) + ((out$CB_out_all[[sample_i]]$HW_CB_eval)[1, , as_right_ind] < 0)) > 0) >= 1)  
      Cao_CB_out01[as_right_ind] = (sum((((out$CB_out_all[[sample_i]]$CaoCB_Nair[[as_right_ind]])[, 1] > parameters$tau) + ((out$CB_out_all[[sample_i]]$CaoCB_Nair[[as_right_ind]])[, 2] < 0)) > 0) >= 1)  
      Cao_CB2_out01[as_right_ind] = (sum((((out$CB_out_all[[sample_i]]$CaoCB2_Nair[[as_right_ind]])[, 1] > parameters$tau) + ((out$CB_out_all[[sample_i]]$CaoCB2_Nair[[as_right_ind]])[, 2] < 0)) > 0) >= 1)  
      Bs_CB_out01[as_right_ind] = (sum((((out$CB_out_all[[sample_i]]$GeoCB_Nair_Bs[[as_right_ind]])[, 1] > parameters$tau) + ((out$CB_out_all[[sample_i]]$GeoCB_Nair_Bs[[as_right_ind]])[, 2] < 0)) > 0) >= 1)  
      BEc_CB_out01[as_right_ind] = (sum((((out$CB_out_all[[sample_i]]$GeoCB_Nair_noBs[[as_right_ind]])[, 1] > parameters$tau) + ((out$CB_out_all[[sample_i]]$GeoCB_Nair_noBs[[as_right_ind]])[, 2] < 0)) > 0) >= 1)  
      naivet_CB_out01[as_right_ind] = (sum((((out$CB_out_all[[sample_i]]$GeoCB_Nair_noBs[[as_right_ind]])[, 3] > parameters$tau) + ((out$CB_out_all[[sample_i]]$GeoCB_Nair_noBs[[as_right_ind]])[, 4] < 0)) > 0) >= 1)  
      MFD_CB_out01[as_right_ind] = (sum((((out$CB_out_all[[sample_i]]$MFD_CB_Nair[[as_right_ind]])[2, ] > parameters$tau) + ((out$CB_out_all[[sample_i]]$MFD_CB_Nair[[as_right_ind]])[1, ] < 0)) > 0) >= 1)  
      CB_out01_mat[parameters$nrep_sub * (split_ind - 1) + sample_i, , as_right_ind] = c(EL_CB_Nair_out01[as_right_ind], EP_CB_Nair_out01[as_right_ind], HW_CB_out01[as_right_ind], Cao_CB_out01[as_right_ind], Cao_CB2_out01[as_right_ind], Bs_CB_out01[as_right_ind], BEc_CB_out01[as_right_ind], naivet_CB_out01[as_right_ind], MFD_CB_out01[as_right_ind]) 
      EL_CB_Nair_length[as_right_ind] = mean((out$CB_out_all[[sample_i]]$EL_CB_Nair[[as_right_ind]])[2, ] - (out$CB_out_all[[sample_i]]$EL_CB_Nair[[as_right_ind]])[1, ])
      EP_CB_Nair_length[as_right_ind] = mean((out$CB_out_all[[sample_i]]$EP_CB_Nair[[as_right_ind]])[2, ] - (out$CB_out_all[[sample_i]]$EP_CB_Nair[[as_right_ind]])[1, ])
      HW_CB_length[as_right_ind] = mean((out$CB_out_all[[sample_i]]$HW_CB_eval)[2, , as_right_ind] - (out$CB_out_all[[sample_i]]$HW_CB_eval)[1, , as_right_ind])
      Cao_CB_length[as_right_ind] = mean((out$CB_out_all[[sample_i]]$CaoCB_Nair[[as_right_ind]])[, 1] - (out$CB_out_all[[sample_i]]$CaoCB_Nair[[as_right_ind]])[, 2])
      Cao_CB2_length[as_right_ind] = mean((out$CB_out_all[[sample_i]]$CaoCB2_Nair[[as_right_ind]])[, 1] - (out$CB_out_all[[sample_i]]$CaoCB2_Nair[[as_right_ind]])[, 2])
      Bs_CB_length[as_right_ind] = mean((out$CB_out_all[[sample_i]]$GeoCB_Nair_Bs[[as_right_ind]])[, 1] - (out$CB_out_all[[sample_i]]$GeoCB_Nair_Bs[[as_right_ind]])[, 2])
      BEc_CB_length[as_right_ind] = mean((out$CB_out_all[[sample_i]]$GeoCB_Nair_noBs[[as_right_ind]])[, 1] - (out$CB_out_all[[sample_i]]$GeoCB_Nair_noBs[[as_right_ind]])[, 2])
      naivet_CB_length[as_right_ind] = mean((out$CB_out_all[[sample_i]]$GeoCB_Nair_noBs[[as_right_ind]])[, 3] - (out$CB_out_all[[sample_i]]$GeoCB_Nair_noBs[[as_right_ind]])[, 4])
      MFD_CB_length[as_right_ind] = mean((out$CB_out_all[[sample_i]]$MFD_CB_Nair[[as_right_ind]])[2, ] - (out$CB_out_all[[sample_i]]$MFD_CB_Nair[[as_right_ind]])[1, ])
      CB_lengths_mean_mat[parameters$nrep_sub * (split_ind - 1) + sample_i, , as_right_ind] = c(EL_CB_Nair_length[as_right_ind], EP_CB_Nair_length[as_right_ind], HW_CB_length[as_right_ind], Cao_CB_length[as_right_ind], Cao_CB2_length[as_right_ind], Bs_CB_length[as_right_ind], BEc_CB_length[as_right_ind], naivet_CB_length[as_right_ind], MFD_CB_length[as_right_ind])     
      EL_Nair_mono[as_right_ind] = (sum(diff((out$CB_out_all[[sample_i]]$EL_CB_Nair[[as_right_ind]])[2, ] ) > 0) == 0) + (sum(diff((out$CB_out_all[[sample_i]]$EL_CB_Nair[[as_right_ind]])[1, ] ) > 0) == 0)
      EP_Nair_mono[as_right_ind] = (sum(diff((out$CB_out_all[[sample_i]]$EP_CB_Nair[[as_right_ind]])[2, ] ) > 0) == 0) + (sum(diff((out$CB_out_all[[sample_i]]$EP_CB_Nair[[as_right_ind]])[1, ] ) > 0) == 0)
      HW_mono[as_right_ind] = (sum(diff((out$CB_out_all[[sample_i]]$HW_CB_eval)[2, , as_right_ind]) > 0) == 0) + (sum(diff((out$CB_out_all[[sample_i]]$HW_CB_eval)[1, , as_right_ind]) > 0) == 0)
      Cao_mono[as_right_ind] = (sum(diff((out$CB_out_all[[sample_i]]$Cao_CB_Nair[[as_right_ind]])[, 1]) > 0) == 0) + (sum(diff((out$CB_out_all[[sample_i]]$Cao_CB_Nair[[as_right_ind]])[, 2]) > 0) == 0)
      Cao2_mono[as_right_ind] = (sum(diff((out$CB_out_all[[sample_i]]$Cao_CB2_Nair[[as_right_ind]])[, 1]) > 0) == 0) + (sum(diff((out$CB_out_all[[sample_i]]$Cao_CB2_Nair[[as_right_ind]])[, 2]) > 0) == 0)
      Bs_mono[as_right_ind] = (sum(diff((out$CB_out_all[[sample_i]]$GeoCB_Nair_Bs[[as_right_ind]])[, 1]) > 0) == 0) + (sum(diff((out$CB_out_all[[sample_i]]$GeoCB_Nair_Bs[[as_right_ind]])[, 2]) > 0) == 0)
      BEc_mono[as_right_ind] = (sum(diff((out$CB_out_all[[sample_i]]$GeoCB_Nair_noBs[[as_right_ind]])[, 1]) > 0) == 0) + (sum(diff((out$CB_out_all[[sample_i]]$GeoCB_Nair_noBs[[as_right_ind]])[, 2]) > 0) == 0)
      naivet_mono[as_right_ind] = (sum(diff((out$CB_out_all[[sample_i]]$GeoCB_Nair_noBs[[as_right_ind]])[, 3]) > 0) == 0) + (sum(diff((out$CB_out_all[[sample_i]]$GeoCB_Nair_noBs[[as_right_ind]])[, 4]) > 0) == 0)
      MFD_mono[as_right_ind] = (sum(diff((out$CB_out_all[[sample_i]]$MFD_CB_Nair[[as_right_ind]])[2, ]) > 0) == 0) + (sum(diff((out$CB_out_all[[sample_i]]$MFD_CB_Nair[[as_right_ind]])[1, ]) > 0) == 0)
      CB_mono_mat[parameters$nrep_sub * (split_ind - 1) + sample_i, , as_right_ind] = c(EL_Nair_mono[as_right_ind], EP_Nair_mono[as_right_ind], HW_mono[as_right_ind], Cao_mono[as_right_ind], Cao2_mono[as_right_ind], Bs_mono[as_right_ind], BEc_mono[as_right_ind], naivet_mono[as_right_ind], MFD_mono[as_right_ind]) 
    }
  }
  split_ind = split_ind + 1
}  
# computing the empirical coverage rates:
cov_result = matrix(apply(cov_mat, 2, mean), nrow = right_endpt_len)  
colnames(cov_result) = c('EL_Nair', 'EP_Nair', 'HW', 'Cao', 'Cao2', 'Bs', 'BEc', 'naivet', 'MFD')
cov_result
# computing the average widths:
t(apply(CB_lengths_mean_mat, c(2, 3), mean))
# computing the range-violation rate:
t(apply(CB_out01_mat, c(2, 3), sum)) 
# computing the average number of confidence band boundaries that satisfy monotonicity (monotone: 2; not: < 2)
t(apply(CB_mono_mat, c(2, 3), mean)) 



