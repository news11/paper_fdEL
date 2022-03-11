# ------------------------------------------------------------------------------------------------------------------ 
# Section 3. Simulation Study - 3.2. Performance of ANOVA Tests
# ------------------------------------------------------------------------------------------------------------------

# ------------ load the packages -------------- 
rm(list=ls())
if (!require(parallel)) {
  install.packages("parallel")
  require(parallel) # this package is needed for the \code{makeCluster} and \code{stopCluster} functions below for cluster computing
}
if (!require(foreach)) {
  install.packages("foreach")
  require(foreach) # we use version 1.5.2 of this package, which is needed for the \code{foreach} function below for cluster computing
}
if (!require(doParallel)) {
  install.packages("doParallel")
  require(doParallel) # we use version 1.0.17 of this package, which is needed for the \code{registerDoParallel} function below for cluster computing
}
if (!require(ESGtoolkit)) {
  install.packages("E:\\R_package\\ESGtoolkit_0.1.tar.gz", repos = NULL, type="source") # change the link to the folder where the downloaded package is; remember to install its dependencies as well
  require(ESGtoolkit) # we use version 0.1 of this package, which is needed for the \code{simdiff} function below for generating \eqn{X(t)}
}
if (!require(devtools)) {
  install.packages("devtools") # we use version 2.4.3 of this package, which is needed for installing the fdEL package below
  require(devtools)
}
if (!require(fdEL)) {
  install_github("news11/fdEL")
  require(fdEL) # we use version 0.0.0.9000 of this package, which is needed for implementing the EL procedures, the Cao and Cao2 confidence bands, and the Geo and Fmax tests
} 
if (!require(fdANOVA)) {
  install.packages("fdANOVA")
  require(fdANOVA) # we use version 0.1.2 of this package, which is needed for implementing the test TRP in the function fanova.tests
}

# ------------ user input -------------- #
dir_path2 = "/home/hwchang/NA/rdata" # where the resulting files can be saved
parameters = list()
# Results in Table 2, first row correspond to setting n_subject = c(70, 100, 130), theta3sq_div_theta2 =  c(64/24, 64/24, 64/24), theta2 = c(6, 6, 6), theta1_div_theta2 = c(-0.4, -0.4, -0.4), beta1_div_beta2 = c(10, 10, 10), beta2 = c(16, 1.1, 0.09) below
# Results in Table 2, second row correspond to setting n_subject = c(130, 100, 70), theta3sq_div_theta2 =  c(64/24, 64/24, 64/24), theta2 = c(6, 6, 6), theta1_div_theta2 = c(-0.4, -0.4, -0.4), beta1_div_beta2 = c(10, 10, 10), beta2 = c(16, 1.1, 0.09) below
# Results in Table 2, third row correspond to setting n_subject = c(70, 100, 130), theta3sq_div_theta2 =  c(64/24, 64/24 * 1.024, 64/24 * 1.07), theta2 = c(2, 0.5, 0.225), theta1_div_theta2 = c(-0.4, -0.3916, -0.415), beta1_div_beta2 = c(10, 10, 10), beta2 = c(10, 2, 0.125) below
# Results in Table 2, fourth row correspond to setting n_subject = c(130, 100, 70), theta3sq_div_theta2 =  c(64/24, 64/24 * 1.024, 64/24 * 1.07), theta2 = c(2, 0.5, 0.225), theta1_div_theta2 = c(-0.4, -0.3916, -0.415), beta1_div_beta2 = c(10, 10, 10), beta2 = c(10, 2, 0.125) below
# Results in Table 2, fifth row correspond to setting n_subject = c(70, 100, 130), theta3sq_div_theta2 =  c(64/24, 64/24 * 1.039, 64/24 * 1.045), theta2 = c(2, 0.5, 0.07), theta1_div_theta2 = c(-0.4, -0.395, -0.38), beta1_div_beta2 = c(10, 10, 10), beta2 = c(10, 2, 0.45) below
# Results in Table 2, sixth row correspond to setting n_subject = c(130, 100, 70), theta3sq_div_theta2 =  c(64/24, 64/24 * 1.039, 64/24 * 1.045), theta2 = c(2, 0.5, 0.07), theta1_div_theta2 = c(-0.4, -0.395, -0.38), beta1_div_beta2 = c(10, 10, 10), beta2 = c(10, 2, 0.45) below
parameters$no_cores <- 50 # the number of cores for parallel computing, preferably < parallel::detectCores(); changing it can still reproduce the results in case 2 of table 1, as long as parameters$no_cores * parameters$nrep_sub = 1000
parameters$nrep_sub = 20 # number of datasets to be generated per parallel task; changing it can still reproduce the results in case 2 of table 1, as long as parameters$no_cores * parameters$nrep_sub = 1000
parameters$n_subject = c(130, 100, 70) # a vector whose j-th element is the number of subjects in the j-th group
parameters$nboot = 1000 # number of bootstrap samples 
parameters$alpha_vec = 0.05 # a number representing the significance level of interest
parameters$tau = 1 # [0, tau): the time domain of the \eqn{X(t)} process
parameters$by_agrid = 2 # a number with values >= 1, representing the mesh of \eqn{{\bf G}_n} (defined in Section 2.2, from the smallest to the largest)
parameters$gridsol = 0.001 # the mesh of the domain of \eqn{X(t)} in numerical simulation
# the following are to determine the parameters \eqn{\theta_1, \theta_2, \theta_3} for the Ornstein--Uhlenbeck process \eqn{\Omega(t)} in Section 3.2 (we omit the subscript \eqn{j} for now), where \eqn{d\Omega(t) = (\theta_1 - \theta_2 \Omega(t))dt + \theta_3 dW(t)}, and \eqn{W(t)} is a standard Brownian motion 
parameters$theta3sq_div_theta2 =  c(64/24, 64/24 * 1.024, 64/24 * 1.07) # a vector whose j-th element is \eqn{(\theta3 ^ 2) / \theta2} for the j-th group
parameters$theta2 = c(2, 0.5, 0.225) # a vector whose j-th element is \eqn{\theta_2} for the j-th group
parameters$theta1_div_theta2 = c(-0.4, -0.3916, -0.415) # a vector whose j-th element is \eqn{\theta_1 / \theta_2} for the j-th group
parameters$theta1 = parameters$theta1_div_theta2 * parameters$theta2 # a vector whose j-th element is \eqn{\theta_1} for the j-th group
parameters$theta3 = sqrt(parameters$theta3sq_div_theta2 * parameters$theta2) # a vector whose j-th element is \eqn{\theta_3} for the j-th group
# the following are to determine the parameters \eqn{\beta_1, \beta_2} for the beta random variable \eqn{\Sigma_j} in Section 3.2
parameters$beta1_div_beta2 = c(10, 10, 10)  # a vector whose j-th element is \eqn{\beta_1 / \beta_2} for the j-th group
parameters$beta2 = c(10, 2, 0.125) # a vector whose j-th element is \eqn{\beta_2} for the j-th group
parameters$beta1 = parameters$beta1_div_beta2 * parameters$beta2 # a vector whose j-th element is \eqn{\beta_1} for the j-th group
parameters$right_endpt = c(0.05, 0.025, 0) # the \eqn{{\bf Z}} set in Supplement Section 5.2

# ------------ define functions to be used -------------- #
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
powerfn = function(parameters, n_subject, split) {
  # Simulate data from the example in Section 3.2 and computes the empirical rejection rates of the functional ANOVA tests based on EL and other methods
  # Args:
  #   parameters: a list of parameters from the ``user input'' above
  #   n_subject: a vector whose j-th element is the number of subjects in the j-th group
  #   split: number of parallel tasks
  # Returns: 
  #   a list containing the following:
  #   rej_rates: the empirical rejection rates of the functional ANOVA tests based on EL and other methods
  t_ncol = parameters$tau / parameters$gridsol + 1
  right_endpt_len = length(parameters$right_endpt)
  # define output
  anova_onefactor_rej_all = array(0, c(parameters$nrep_sub, length(parameters$alpha_vec), 3))
  suptest_rej_all = matrix(0, nrow = parameters$nrep_sub, ncol = length(parameters$alpha_vec))
  suptest_EP_rej_all = matrix(0, nrow = parameters$nrep_sub, ncol = length(parameters$alpha_vec))
  # for loop generating \code{parameters$nrep_sub} datasets and computing the test statistics accordingly
  for (nrepi in 1:parameters$nrep_sub) {
    set.seed(parameters$nrep_sub * (split - 1) + nrepi, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rounding") # for R versions < 3.6.0, may need to remove the sample.kind argument to avoid error messages
    # simulate \eqn{X_j(t)} in the example in Section 3.2
    beta_to_time = list()
    Xt = lapply(1:length(n_subject), FUN = function (j) {
      sim.OU <- ESGtoolkit::simdiff(n = n_subject[j], horizon = (t_ncol - 1), frequency = "annual", model = "OU", x0 = 0, theta1 = parameters$theta1[j], theta2 = parameters$theta2[j], theta3 = parameters$theta3[j])
      beta_to_time[[j]] <<- rbeta(n_subject[j], parameters$beta1[j], parameters$beta2[j])
      floor_n(300 * pmax(t(sim.OU), 0))
    })  
    # compute \eqn{L_j(a)} in the example in Section 3.2
    as = seq(0, max(unlist(Xt)), by = parameters$by_agrid)
    La = lapply(1:length(n_subject), FUN = function (j) {
      t(sapply(1:n_subject[j], FUN = function (i) {
        observed_Xt = (Xt[[j]])[i, -1][!is.na((Xt[[j]])[i, -1])] 
        observed_Xt_len = length(observed_Xt)
        grid_width_i = rep(1 / observed_Xt_len, observed_Xt_len)  
        fdEL::mu_a_hat(observed_Xt, grid_widths = grid_width_i, as) * (beta_to_time[[j]])[i]
      }))
    })
    mu_a0 = list()  
    right_ind_keep_as = sapply(1:length(n_subject), FUN = function (j) {
      mu_a0[[j]] <<- apply(La[[j]], 2, mean)
      sapply(parameters$right_endpt, function(i) {
        if (sum(mu_a0[[j]] / parameters$tau > i) != 0) {
          max(which(mu_a0[[j]] / parameters$tau > i))
        } else {1}
      })  
    }) 
    right_ind_keep_as0 = apply(right_ind_keep_as, 1, min)
    # computing the EL and related functional ANOVA tests
    test_out = fdEL::elfanova(Ta = La, as = as, n_boot = parameters$nboot, as_right_ind_keep_as0 = right_ind_keep_as0)
    suptest_rej_all[nrepi, ] = as.numeric(test_out$out_sup_pval < parameters$alpha_vec)
    suptest_EP_rej_all[nrepi, ] = as.numeric(test_out$out_sup_EP_pval < parameters$alpha_vec)
    La_mat = matrix(unlist(lapply(1:length(n_subject), FUN = function (j) {  # for fda.usc's fdata object construction
      c(t(La[[j]]))
    })), byrow = T, ncol = length(as))
    anova_onefactor_group = as.factor(rep(1:length(n_subject), times = n_subject))
    # computing the GPF and Fmax tests mentioned in Section 3
    fanova_FmaxbGPF_Nair <- fdEL::fanovatests_sel(x = t(La_mat), group.label = anova_onefactor_group, n_boot = parameters$nboot, as_right_ind_keep_as0 = right_ind_keep_as0, parameters = list(dir_path2 = getwd(), no_cores = 1, useseed = FALSE, seed = NULL)) 
    # computing the TRP test mentioned in Section 3
    fanova <- fdANOVA::fanova.tests(x = t(La_mat), group.label = anova_onefactor_group, test = c("TRP"), params = list(paramTRP = list(B.TRP = parameters$nboot))) 
    test_out$out_anova_onefactor_pval = c(fanova_FmaxbGPF_Nair, fanova$TRP$pvalues.WTPS)
    anova_onefactor_rej_all[nrepi, 1, ] = as.numeric(test_out$out_anova_onefactor_pval < parameters$alpha_vec[1])
  }  
  return(list(rej_rates = c(mean(suptest_rej_all[, 1]), mean(suptest_EP_rej_all[, 1]), apply(anova_onefactor_rej_all[, 1, ], 2, mean))
  ))   
}  

# ------------ calculate rejection rates by utilizing R parallel computing -------------- #
# initiate cluster
cl = parallel::makeCluster(parameters$no_cores)
doParallel::registerDoParallel(cl)
# compute the results:
foreach::foreach(split = 1:parameters$no_cores,  
  .combine = c,
  .packages = c("parallel", "foreach", "doParallel", "ESGtoolkit", "devtools", "fdEL", "fdANOVA")
) %dopar% {
  out = powerfn(parameters, n_subject = parameters$n_subject, split = split)
  setwd(dir_path2)
# 20220308 parameters$n_floor = 1, parameters$SDEmodel = "OU": save(out,  file = paste("n_", paste(parameters$n_subject, collapse = "_"), "_beta2_", paste(parameters$beta2, collapse = "_"), '_byagrid_', parameters$by_agrid, '_nfloor_', parameters$n_floor, "_split_", split, '_nrep_', parameters$nrep, '_endpt_', paste(parameters$right_endpt, collapse = "_"), '_SDEmodel_', parameters$SDEmodel, "_CIR1_", paste(parameters$CIR1, collapse = "_"), "_CIR2_", paste(parameters$CIR2, collapse = "_"), "_CIR3_", paste(parameters$CIR3, collapse = "_"), ".Rdata", sep = ""))
  save(out,  file = paste("n_", paste(parameters$n_subject, collapse = "_"), "_beta2_", paste(parameters$beta2, collapse = "_"), '_byagrid_', parameters$by_agrid, "_split_", split, '_nrep_sub_', parameters$nrep_sub, '_endpt_', paste(parameters$right_endpt, collapse = "_"), "_theta1_", paste(parameters$theta1, collapse = "_"), "_theta2_", paste(parameters$theta2, collapse = "_"), "_theta3_", paste(parameters$theta3, collapse = "_"), ".Rdata", sep = ""))
}  
parallel::stopCluster(cl)

# ------------ collecting results from R parallel computing -------------- #
t_ncol = parameters$tau / parameters$gridsol + 1
Set = 1:parameters$no_cores
rej_mat = matrix(0, nrow = length(Set), ncol = 5)  # 0807 bug correction
setwd(dir_path2)
split_ind = 1
for (split in Set) {
  load(paste("n_", paste(parameters$n_subject, collapse = "_"), "_beta2_", paste(parameters$beta2, collapse = "_"), '_byagrid_', parameters$by_agrid, "_split_", split, '_nrep_sub_', parameters$nrep_sub, '_endpt_', paste(parameters$right_endpt, collapse = "_"), "_theta1_", paste(parameters$theta1, collapse = "_"), "_theta2_", paste(parameters$theta2, collapse = "_"), "_theta3_", paste(parameters$theta3, collapse = "_"), ".Rdata", sep = ""))
  rej_mat[split_ind, ] = out$rej_rates
  split_ind = split_ind + 1
}  
colnames(rej_mat)=c('supEL','supWald', 'GPF','Fmaxb','TRP_WTPS')
# computing the empirical rejection rates:
apply(rej_mat, 2, mean)  
