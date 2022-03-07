# ------------------------------------------------------------------------------------------------------------------ 
# Section 3. Simulation Study - 3.1. Performance of Simultaneous Confidence Bands - first example
# ------------------------------------------------------------------------------------------------------------------

# ------------ load the packages -------------- #
rm(list=ls())
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)  # create personal library
.libPaths(Sys.getenv("R_LIBS_USER"))  # add to the path
if (!require(parallel)) {
  install.packages("parallel")
  require(parallel) # this package is needed for the \code{makeCluster} and \code{stopCluster} functions below for cluster computing
}
if (!require(foreach)) {
  install.packages("foreach")
  require(foreach) # this package is needed for the \code{foreach} function below for cluster computing
}
if (!require(doParallel)) {
  install.packages("doParallel")
  require(doParallel) # this package is needed for the \code{registerDoParallel} function below for cluster computing
}
if (!require(MASS)) {
  install.packages("MASS")
  require(MASS) # this package is needed for the \code{mvrnorm} function below for generating T(a)
}
#if (!require(fdEL)) {
  install.packages("/home/hwchang/NA/sh/fdEL" , repos = NULL, type = "source")
  require(fdEL) # this package is needed for implementing the EL procedures, the Cao and Cao2 confidence bands, and the Geo and Fmax tests
#} 
if (!require(SCBmeanfd)) {
  install.packages("SCBmeanfd")
  require(SCBmeanfd) # this package is needed for implementing the MFD confidence band
}
if (!require(KernSmooth)) {
  install.packages("KernSmooth")
  require(KernSmooth) # this package is needed for the locpoly function inside the function scb.mean00
}
if (!require(devtools)) {
  install.packages("devtools") # this package is needed for installing the fregion package below
  require(devtools)
}
if (!require(fregion)) {
  install_github("hpchoi/fregion", force = TRUE)
  require(fregion) # this package is needed for implementing the Geo confidence band
}

# ------------ user input -------------- #
dir_path = "/home/hwchang/NA/sh" # where the source files are
dir_path2 = "/home/hwchang/NA/rdata" # where the resulting files can be saved
## set parameters for generating data: (later can be subject specific)
parameters = list()
parameters$no_cores <- 50 # the number of cores for parallel computing
parameters$nrep_sub = 20 # number of datasets to be generated per parallel task
parameters$nboot = 1000 # number of bootstrap samples 
parameters$alpha_vec = 0.05 # a number representing the significance level of interest
parameters$n_subject = 100 # number of subjects
parameters$alpha2 = 1 # [0, alpha2]: the domain of the T(a) process 
parameters$gridsol = 0.001 # the mesh of the domain of T(a) in numerical simulation; this needs to < \code{parameters$alpha2 / parameters$n_agrid_eval} (defined below) and < \code{parameters$alpha2 / parameters$n_agrid} (defined below) 
parameters$n_agrid = 25 + 1 # the number of the design time points of the observed functional data; that is, \eqn{{\bf G}_n} defined in Section 2.2 (from the smallest to the largest).
parameters$n_agrid_eval = 100 + 1 # the number of the time points (from the smallest to the largest; e.g., the activity levels of interest) at which we want to evaluate the confidence band. 
# the following are the three parameters for determing the true covariance function of the Gaussian process T(a) in the first example in Section 3.1
parameters$sigma_cs = 0.6 # a number representing \eqn{\mbox{cov}(T(a), T(a))} for \eqn{a \ge 0.25}
parameters$sigma_jump = 1 # the difference between \eqn{\mbox{cov}(T(a), T(a)) (a \ne b)}  for \eqn{a, b < 0.25} and the case for \eqn{a or b \ge 0.25}
parameters$sigma = 14.6 # a number with values 9.6 or 14.6, representing \eqn{\mbox{cov}(T(a), T(a))} for \eqn{a < 0.25} - \code{parameters$sigma_jump}
parameters$mu_t_b = 0 # a number giving the value of the true functional mean at design time points < 0.25

# ------------ define functions to be used -------------- #
mu_t_fn = function(t_vec, mu_t_b = parameters$mu_t_b) {
  # The true mean of the Gaussian process in the first example in Section 3.1
  # Args:
  #   t_vec: a vector of design time points at which the mean function is evaluated 
  #   mu_t_b: a number giving the value of the true functional mean at \code{t_vec}
  # Returns: 
  #   the desired mean function at \code{t_vec} 
  mu_vec = 1:length(t_vec) * 0
  mu_vec[1:length(t_vec)] = mu_t_b
  return(mu_vec)
}
sigma_t_fn = function(t_vec) {
  # The true covariance function of the Gaussian process in the first example in Section 3.1
  # Args:
  #   t_vec: a vector of design time points of interest 
  # Returns: 
  #   a matrix containing the desired covariance function at each pair of \code{t_vec}
  t_vec_len = length(t_vec)
  cs_mat = toeplitz(c(parameters$sigma_cs, rep(0.5, length = t_vec_len - 1)))  
  diag(cs_mat)[t_vec < 0.25] = parameters$sigma
  jump_mat = matrix(0,  t_vec_len,  t_vec_len)
  jump_mat[t_vec < 0.25, t_vec < 0.25] = parameters$sigma_jump
  sigma_mat = cs_mat + jump_mat
  return(sigma_mat)
}

setwd(dir_path) # for sourcing files in the following lines
source('scb.mean00.R', local = TRUE)
source('Geoband.R', local = TRUE)
coveragefn = function(parameters, n_subject, split) { 
  # Simulate data from the first example in Section 3.1 and computes the empirical coverage of EL, EP and NS confidence bands
  # Args:
  #   parameters: a list of parameters from the ``user input'' above
  #   n_subject: the number of subjects 
  #   split: number of parallel tasks
  # Returns: 
  #   a list containing the following:
  #   CB_out_all: a list whose \code{nrepi}-th element contains some of the outputs from the functions for computing the EL, EP, NS, Cao, Cao2, MFDbs, Geo, naivet and MFD confidence bands based on the \code{nrepi}-th simulated dataset
  #   coverages: the empirical coverage of EL, EP, NS, Cao, Cao2, MFDbs, Geo, naivet and MFD confidence bands
  t_ncol = parameters$alpha2 / parameters$gridsol  
  # define output
  EL_CB_Nairrej_all = 1:parameters$nrep_sub * 0
  EP_CB_Nairrej_all = 1:parameters$nrep_sub * 0
  HW_CBrej_all = 1:parameters$nrep_sub * 0
  CB_out_all = list()
  Cao_CB0rej_all = 1:parameters$nrep_sub * 0
  Cao_CBrej_all = 1:parameters$nrep_sub * 0
  Cao_CB2rej_all = 1:parameters$nrep_sub * 0
  Geo_CBrej_all = matrix(0, nrow = parameters$nrep_sub, ncol = 3) 
  MFD_CBrej_all = 1:parameters$nrep_sub * 0
  # for loop generating \code{parameters$nrep_sub} datasets and computing the confidence bands accordingly
  for (nrepi in 1:parameters$nrep_sub) {  
    set.seed(parameters$nrep_sub * (split - 1) + nrepi, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rounding") 
    t_s_as = seq(0, 1, length = parameters$n_agrid)
    t_s_eval = seq(0, 1, length = parameters$n_agrid_eval)
    mu_t_agrid = mu_t_fn(t_s_as)
    mu_t_eval = mu_t_fn(t_s_eval)
    sigma_t_agrid = sigma_t_fn(t_s_as)
    sigma_t_eval = sigma_t_fn(t_s_eval)
    mu_t_eval_new = mu_t_eval
    # compute the true functional mean in the first example in Section 3.1
    mu_t_eval_new[mu_t_eval >= 0] = mu_t_eval[mu_t_eval >= 0] + sqrt(diag(sigma_t_eval)[mu_t_eval >= 0]) / sqrt(2 * pi) - mu_t_eval[mu_t_eval >= 0] * pnorm(-mu_t_eval[mu_t_eval >= 0] / sqrt(diag(sigma_t_eval)[mu_t_eval >= 0]))
    mu_t_eval_new[mu_t_eval < 0] = mu_t_eval[mu_t_eval < 0] + sqrt(diag(sigma_t_eval)[mu_t_eval < 0]) * exp(-(mu_t_eval[mu_t_eval < 0] ^ 2)/ (2 * diag(sigma_t_eval)[mu_t_eval < 0])) / sqrt(2 * pi) - mu_t_eval[mu_t_eval < 0] * pnorm(-mu_t_eval[mu_t_eval < 0] / sqrt(diag(sigma_t_eval)[mu_t_eval < 0]))
    mu_t_eval = mu_t_eval_new 
    # generate data T(a) in the first example in Section 3.1
    Ta = t(sapply(1:n_subject, FUN = function (i) {
      pmax(MASS::mvrnorm(mu = mu_t_agrid, Sigma = sigma_t_agrid), 0)
    }))
    # computing the EL confidence band
    CB_out = fdEL::elfband(Ta = Ta, as = t_s_as, as_eval = t_s_eval, n_boot = parameters$nboot, mu_a = mu_t_eval, alpha = parameters$alpha_vec)
    EL_CB_Nairrej_all[nrepi] = CB_out$EL_CB_Nairrej
    EP_CB_Nairrej_all[nrepi] = CB_out$EP_CB_Nairrej
    HW_CBrej_all[nrepi] = CB_out$HW_CBrej
    # NOT USED: computing the Cao confidence band before projecting the initially smoothed covariance estimates onto the space of non-negative definite matrices    
    an.error.occured <- FALSE
    a.warning.occured <- FALSE
    possibleError_Cao0 <- tryCatch(
      CaoCB0 <- fdEL::Caoband(Y = Ta, as = t_s_as, checkpts = t_s_eval, b_M = parameters$nboot, mean_true = mu_t_eval, alpha = parameters$alpha_vec[1], cov_nnd = 0),
      error = function(e) {an.error.occured <<- TRUE},
      warning = function(w) {a.warning.occured <<- TRUE}
    ) 
    if (sum(an.error.occured + a.warning.occured) > 0) {
      CB_out$CaoCB0 = matrix(-1, ncol = 2, nrow = parameters$n_agrid_eval) 
      Cao_CB0rej_all[nrepi] = 1
    } else {
      CB_out$CaoCB0 = CaoCB0$confidence_bounds
      Cao_CB0rej_all[nrepi] = CaoCB0$cover_or_not
    }
    # computing the Cao confidence band mentioned in Section 3    
    CaoCB <- fdEL::Caoband(Y = Ta, as = t_s_as, checkpts = t_s_eval, b_M = parameters$nboot, mean_true = mu_t_eval, alpha = parameters$alpha_vec[1])
    CB_out$CaoCB = CaoCB$confidence_bounds
    Cao_CBrej_all[nrepi] = CaoCB$cover_or_not
    # computing the Cao2 confidence band mentioned in Section 3    
    CaoCB2 <- fdEL::Caoband(Y = Ta, as = t_s_as, checkpts = t_s_eval, b_M = parameters$nboot, mean_true = mu_t_eval, alpha = parameters$alpha_vec[1], N_m = parameters$n_agrid - 2)
    CB_out$CaoCB2 = CaoCB2$confidence_bounds
    Cao_CB2rej_all[nrepi] = CaoCB2$cover_or_not
    # computing the MFDbs and Geo confidence bands mentioned in Section 3 (\code{btypes = "naive.t"} below is not reported in Section 3 but for comparison purpose only, because it only guarantees pointwise coverage accuracy instead of simultaneous coverage accuracy)
    GeoCB_Bs <- fmean(Y=t(Ta), as = t_s_as, as_eval = t_s_eval, n_boot = parameters$nboot, alpha = parameters$alpha_vec[1], btypes =  c("Bs"), as_right_ind_keep_as0 = parameters$n_agrid, right_endpt = 0)   
    GeoCB_noBs <- fmean(Y=t(Ta), as = t_s_as, as_eval = t_s_eval, n_boot = parameters$nboot, alpha = parameters$alpha_vec[1], btypes = c("BEc", "naive.t"), as_right_ind_keep_as0 = parameters$n_agrid, right_endpt = 0)   
    CB_out$GeoCB = cbind(GeoCB_Bs$confidence_bounds, GeoCB_noBs$confidence_bounds)  # columns are: Bs.u, Bs.l, BEc.u, BEc.l, naivet.u, naivet.l
    Geo_CBrej_all[nrepi,] = c((any(CB_out$GeoCB[, 1] < mu_t_eval) || any(CB_out$GeoCB[, 2] > mu_t_eval)), (any(CB_out$GeoCB[, 3] < mu_t_eval) || any(CB_out$GeoCB[, 4] > mu_t_eval)), (any(CB_out$GeoCB[, 5] < mu_t_eval) || any(CB_out$GeoCB[, 6] > mu_t_eval)))
    # computing the MFD confidence band mentioned in Section 3
    h <- SCBmeanfd::cv.select(t_s_as, Ta, 1)  
    scbTa <- scb.mean00(t_s_as, Ta, bandwidth = h, degree = 1, gridsize = parameters$n_agrid_eval) 
    CB_out$MFD_CB = t(scbTa$normscb)  # 1st row = lower bound; 2nd row = upper bound 
    MFD_CBrej_all[nrepi] = (any(CB_out$MFD_CB[2, ] < mu_t_eval) || any(CB_out$MFD_CB[1, ] > mu_t_eval))
    CB_out_all[[nrepi]] = CB_out
  } # END for
  return(list(CB_out_all = CB_out_all, coverages = c(mean(1 - EL_CB_Nairrej_all), mean(1 - EP_CB_Nairrej_all), mean(1 - HW_CBrej_all), mean(1 - Cao_CBrej_all), mean(1 - Cao_CB2rej_all), apply(1 - Geo_CBrej_all, 2, mean), mean(1 - MFD_CBrej_all))))   
} # END coveragefn
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
  out = coveragefn(parameters, n_subject = parameters$n_subject, split = split)
  setwd(dir_path2)
  save(out,  file = paste("n_", parameters$n_subject, "_split_", split, '_nrep_sub_', parameters$nrep_sub, "_sigma_", parameters$sigma, "_sigma_jump_", parameters$sigma_jump, "_sigma_cs_", parameters$sigma_cs, "_mu_t_b_", parameters$mu_t_b, "_gridsol_", parameters$gridsol, "_n_agrid_", parameters$n_agrid, "_n_agrid_eval_", parameters$n_agrid_eval, ".Rdata", sep = ""))
}  # END coveragefn
parallel::stopCluster(cl)

# ------------ collecting results from R parallel computing -------------- #
t_ncol = parameters$alpha2 / parameters$gridsol
Set = 1:parameters$no_cores  
cov_mat = matrix(0, nrow = length(Set), ncol = 9)
CB_out01_mat = matrix(0, nrow = length(Set) * parameters$nrep_sub, ncol = 9)
CB_lengths_mean_mat = matrix(0, nrow = length(Set) * parameters$nrep_sub, ncol = 9)
setwd(dir_path2)
split_ind = 1
for (split in Set) {
  load(paste("n_", parameters$n_subject, "_split_", split, '_nrep_sub_', parameters$nrep_sub, "_sigma_", parameters$sigma, "_sigma_jump_", parameters$sigma_jump, "_sigma_cs_", parameters$sigma_cs, "_mu_t_b_", parameters$mu_t_b, "_gridsol_", parameters$gridsol, "_n_agrid_", parameters$n_agrid, "_n_agrid_eval_", parameters$n_agrid_eval, ".Rdata", sep = ""))
  cov_mat[split_ind, ] = out$coverages
  for (sample_i in 1:parameters$nrep_sub) {
    EL_CB_Nair_out01 = (sum((((out$CB_out_all[[sample_i]]$EL_CB_Nair[[1]])[2, ] > Inf) + ((out$CB_out_all[[sample_i]]$EL_CB_Nair[[1]])[1, ] < 0)) > 0) >= 1)  
    EP_CB_Nair_out01 = (sum((((out$CB_out_all[[sample_i]]$EP_CB_Nair[[1]])[2, ] > Inf) + ((out$CB_out_all[[sample_i]]$EP_CB_Nair[[1]])[1, ] < 0)) > 0) >= 1)  
    HW_CB_out01 = (sum((((out$CB_out_all[[sample_i]]$HW_CB)[2, ] > Inf) + ((out$CB_out_all[[sample_i]]$HW_CB)[1, ] < 0)) > 0) >= 1)  
    Cao_CB_out01 = (sum((((out$CB_out_all[[sample_i]]$CaoCB)[, 1] > Inf) + ((out$CB_out_all[[sample_i]]$CaoCB)[, 2] < 0)) > 0) >= 1)  
    Cao_CB2_out01 = (sum((((out$CB_out_all[[sample_i]]$CaoCB2)[, 1] > Inf) + ((out$CB_out_all[[sample_i]]$CaoCB2)[, 2] < 0)) > 0) >= 1)  
    Bs_CB_out01 = (sum((((out$CB_out_all[[sample_i]]$GeoCB)[, 1] > Inf) + ((out$CB_out_all[[sample_i]]$GeoCB)[, 2] < 0)) > 0) >= 1)  
    BEc_CB_out01 = (sum((((out$CB_out_all[[sample_i]]$GeoCB)[, 3] > Inf) + ((out$CB_out_all[[sample_i]]$GeoCB)[, 4] < 0)) > 0) >= 1)  
    naivet_CB_out01 = (sum((((out$CB_out_all[[sample_i]]$GeoCB)[, 5] > Inf) + ((out$CB_out_all[[sample_i]]$GeoCB)[, 6] < 0)) > 0) >= 1)  
    MFD_CB_out01 = (sum((((out$CB_out_all[[sample_i]]$MFD_CB)[2, ] > Inf) + ((out$CB_out_all[[sample_i]]$MFD_CB)[1, ] < 0)) > 0) >= 1)  
    CB_out01_mat[parameters$nrep_sub * (split_ind - 1) + sample_i, ] = c(EL_CB_Nair_out01, EP_CB_Nair_out01, HW_CB_out01, Cao_CB_out01, Cao_CB2_out01, Bs_CB_out01, BEc_CB_out01, naivet_CB_out01, MFD_CB_out01)
    EL_CB_Nair_length = (out$CB_out_all[[sample_i]]$EL_CB_Nair[[1]])[2, ] - (out$CB_out_all[[sample_i]]$EL_CB_Nair[[1]])[1, ]
    EP_CB_Nair_length = (out$CB_out_all[[sample_i]]$EP_CB_Nair[[1]])[2, ] - (out$CB_out_all[[sample_i]]$EP_CB_Nair[[1]])[1, ]
    HW_CB_length = (out$CB_out_all[[sample_i]]$HW_CB)[2, ] - (out$CB_out_all[[sample_i]]$HW_CB)[1, ]
    Cao_CB_length = (out$CB_out_all[[sample_i]]$CaoCB)[, 1] - (out$CB_out_all[[sample_i]]$CaoCB)[, 2]
    Cao_CB2_length = (out$CB_out_all[[sample_i]]$CaoCB2)[, 1] - (out$CB_out_all[[sample_i]]$CaoCB2)[, 2]
    Bs_CB_length = (out$CB_out_all[[sample_i]]$GeoCB)[, 1] - (out$CB_out_all[[sample_i]]$GeoCB)[, 2]
    BEc_CB_length = (out$CB_out_all[[sample_i]]$GeoCB)[, 3] - (out$CB_out_all[[sample_i]]$GeoCB)[, 4]
    naivet_CB_length = (out$CB_out_all[[sample_i]]$GeoCB)[, 5] - (out$CB_out_all[[sample_i]]$GeoCB)[, 6]
    MFD_CB_length = (out$CB_out_all[[sample_i]]$MFD_CB)[2, ] - (out$CB_out_all[[sample_i]]$MFD_CB)[1, ]
    CB_lengths_mean_mat[parameters$nrep_sub * (split_ind - 1) + sample_i, ] = c(mean(EL_CB_Nair_length), mean(EP_CB_Nair_length), mean(HW_CB_length), mean(Cao_CB_length), mean(Cao_CB2_length), mean(Bs_CB_length), mean(BEc_CB_length), mean(naivet_CB_length), mean(MFD_CB_length))
  }
  split_ind = split_ind + 1
} 
# computing the empirical coverage rates:
cov_result = apply(cov_mat, 2, mean)  
names(cov_result) = c('EL_Nair', 'EP_Nair', 'HW', 'Cao', 'Cao2', 'Bs', 'BEc', 'naivet', 'MFD')
cov_result 
# computing the average widths:
apply(CB_lengths_mean_mat, 2, mean)
# computing the range-violation rate:
apply(CB_out01_mat, 2, sum) 





