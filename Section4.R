# ------------------------------------------------------------------------------------------------------------------ 
# Section 4. Data Example
# ------------------------------------------------------------------------------------------------------------------
# ------------ load the packages -------------- #
rm(list=ls())
if (!require(fdEL)) {
  install.packages("E:\\R_package\\fdEL", repos = NULL, type = "source")
  require(fdEL) # this package is needed for implementing the EL procedures, the Cao and Cao2 confidence bands, and the Geo and Fmax tests
} 
if (!require(fdANOVA)) {
  install.packages("fdANOVA")
  require(fdANOVA) # this package is needed for implementing the test TRP in the function fanova.tests
}
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
  install_github("hpchoi/fregion")
  require(fregion) # this package is needed for implementing the Geo confidence band
}
if (!require(plotrix)) {
  install.packages("plotrix")
  require(plotrix) # this package is needed for plotted the right panel of Figure 4
}

# ------------ user input -------------- #

parameters = list()
parameters$by_agrid = 1 # the mesh of the grid of activity levels over which the occupation time is evaluated.
parameters$tau = 1
parameters$right_endpt = c(0.05, 0.025, 0) # the \mathbb{Z} set in Supplement Section 6
parameters$nboot = 1000 # number of bootstrap samples
parameters$alpha_vec = 0.05 # significance level
parameters$dir_path_source = "E:\\R_program\\paper_fnl_codes_20220111\\source" # where the source files are
parameters$dir_path2 = "E:\\R_program\\fnl_mil_20210607" # where the resulting files can be saved

# ------------ building occupation time curve from raw activity data Xt -------------- #

as = seq(0, 499, by = parameters$by_agrid) # the grid of activity levels over which the occupation time is evaluated.
n_agrid = length(as)
rawlens = list()
n_group = length(Xt)
# result 1 in Section 4---the sample size of each group:
(n_subject = sapply(1:n_group, FUN = function (j) {
  length(Xt[[j]])
}))
# the occupation time curve:
Ta = lapply(1:n_group, FUN = function (j) {
  observed_Xt_lens_tmp = 1:n_subject[j] * 0
  Ta_j = t(sapply(1:n_subject[j], FUN = function (i) {
    observed_Xt = (Xt[[j]])[[i]][!is.na((Xt[[j]])[[i]])]
    observed_Xt_len = length(observed_Xt)
    observed_Xt_lens_tmp[i] <<- observed_Xt_len
    grid_width_i = rep(1 / observed_Xt_len, observed_Xt_len)  
    mu_a_hat(observed_Xt, grid_widths = grid_width_i, as)
  }))
  rawlens[[j]] <<- observed_Xt_lens_tmp
  return(Ta_j)
}) # Ta is a list of four elements corresponding to the four groups: veterans aged 75-and-older, non-veterans aged 75-and-older, veterans aged 65--74, and non-veterans aged aged 65--74. Each element contains a matrix of occupation time data, with each row containing the occupation time curve of each subject and each column corresponding to each activity level on the grid over which the occupation time is evaluated. rawlens is a list of four elements corresponding to the four groups: veterans aged 75-and-older, non-veterans aged 75-and-older, veterans aged 65--74, and non-veterans aged aged 65--74. Each element contains a vector recording the length of non-NA raw activity data from each subject.
# result 2 in Section 4---number of incomplete curves in each group:
(n_miss_agemilgroups = sapply(1:n_group, FUN = function (j) {
  sum(rawlens[[j]] != 10080)
}))
# result 3 in Section 4---% of missing readings in groups 2 & 4:
range(c(1 - (rawlens[[2]])[rawlens[[2]] != 10080] / 10080, 1 - (rawlens[[4]])[rawlens[[4]] != 10080] / 10080))
# for later R function input: 
mu_a0 = list()  
right_ind_keep_as_mat = sapply(1:length(n_subject), FUN = function (j) {
  mu_a0[[j]] <<- apply(Ta[[j]], 2, mean)
  sapply(parameters$right_endpt, function(i) {
    if (sum(mu_a0[[j]] / parameters$tau > i) != 0) {
      max(which(mu_a0[[j]] / parameters$tau > i))
    } else {1}
  }) 
}) # j-th column is for j-th sample; i-th row is for right_endpt[i]
right_ind_keep_as0 = apply(right_ind_keep_as_mat, 1, min) # a vector; the i-th element is the index for `as' corresponding to \hat{r}(z), z being the i-th element in parameters$right_endpt (i.e., \mathbb{Z})
# another format for occupation time curves for later input related to the fdANOVA package: 
Ta_mat = matrix(unlist(lapply(1:length(n_subject), FUN = function (j) { 
  c(t(Ta[[j]]))
})), byrow = T, ncol = n_agrid)
anova_onefactor_group = as.factor(rep(1:length(n_subject), times = n_subject))

# ------------ functional ANOVA -------------- #

# result 4 in Section 4---numbers in the `all groups' row of Table 3:
set.seed(1008, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
test_out = elfanova(Ta = Ta, as = as, n_boot = parameters$nboot, as_right_ind_keep_as0 = right_ind_keep_as0)
c(test_out$out_sup_pval, test_out$out_sup_EP_pval) # p-values for comparing all 4 groups using the EL and the Wald-type tests, respectively; = 0 means `p < 1 / parameters$nboot'
(fanova_FmaxbGPF_Nair <- fanovatests_sel(x = t(Ta_mat), group.label = anova_onefactor_group, n_boot = parameters$nboot, as_right_ind_keep_as0 = right_ind_keep_as0)) # p-values for comparing all 4 groups using the the GPF and Fmax tests, respectively, which are modified from the GPF and Fmaxb (bootstrapped version of Fmax test) options of the function fanova.tests in the package fdANOVA to implement the Uno's selection approach in Supp Sec 4.2; == 0 means `p < 1 / parameters$nboot'
set.seed(1009, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
(fanova <- fanova.tests(x = t(Ta_mat), group.label = anova_onefactor_group, test = c("TRP"), params = list(paramTRP = list(B.TRP = parameters$nboot)))) 
fanova$TRP$pvalues.WTPS # p-value for comparing all 4 groups using the WTPS test among the TRP tests in the function fanova.tests in the package fdANOVA; p-value == 0 means `p < 1 / parameters$nboot'
# result 5 in Section 4---numbers from the second to the last rows of Table 3:
group_sub = c(3, 4) # change this to c(3, 4), c(1, 3), c(2, 4) for different pairwise comparisons among the groups
set.seed(1008+sum(group_sub), kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
test_out_sub = elfanova(Ta = Ta[group_sub], as = as, n_boot = parameters$nboot, as_right_ind_keep_as0 = right_ind_keep_as0)
c(test_out_sub$out_sup_pval, test_out_sub$out_sup_EP_pval) # p-values for comparing groups in the group_sub vector using the EL and the Wald-type tests, respectively; = 0 means `p < 1 / parameters$nboot'
group_sub = c(1, 2) # change this to c(3, 4), c(1, 3), c(2, 4) for different pairwise comparisons among the groups
Ta_mat_sub = rbind(Ta_mat[anova_onefactor_group == group_sub[1], ], Ta_mat[anova_onefactor_group == group_sub[2], ])
anova_onefactor_group_sub = c(anova_onefactor_group[anova_onefactor_group %in% group_sub[1]], anova_onefactor_group[anova_onefactor_group %in% group_sub[2]])
(fanova1_FmaxbGPF_Nair_sub <- fanovatests_sel(x = t(Ta_mat_sub), group.label = anova_onefactor_group_sub, n_boot = parameters$nboot, as_right_ind_keep_as0 = right_ind_keep_as0, parameters = list(dir_path2 = getwd(), no_cores = 40, useseed = TRUE, seed = 1007))) # p-values for comparing groups in the group_sub vector using the the GPF and Fmax tests, respectively, which are modified from the GPF and Fmaxb (bootstrapped version of Fmax test) options of the function fanova.tests in the package fdANOVA to implement the Uno's selection approach in Supp Sec 4.2; == 0 means `p < 1 / parameters$nboot'
group_sub = c(1, 2) # change this to c(3, 4), c(1, 3), c(2, 4) for different pairwise comparisons among the groups
Ta_mat_sub = rbind(Ta_mat[anova_onefactor_group == group_sub[1], ], Ta_mat[anova_onefactor_group == group_sub[2], ])
anova_onefactor_group_sub = c(anova_onefactor_group[anova_onefactor_group %in% group_sub[1]], anova_onefactor_group[anova_onefactor_group %in% group_sub[2]])
set.seed(1008+sum(group_sub), kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
fanova1 <- fanova.tests(x = t(Ta_mat_sub), group.label = anova_onefactor_group_sub, test = c("TRP"), params = list(paramTRP = list(B.TRP = parameters$nboot))) 
fanova1$TRP$pvalues.WTPS # p-value for comparing groups in the group_sub vector using the WTPS test among the TRP tests in the function fanova.tests in the package fdANOVA; p-value == 0 means `p < 1 / parameters$nboot'

# ------------ confidence band for functional mean -------------- #

setwd(parameters$dir_path_source) # for sourcing files in the following lines
source('scb.mean00.R', local = TRUE)
source('Geoband.R', local = TRUE)
as_eval = as
right_endpt = parameters$right_endpt[1] # here \code{right_endpt} is a number (see \eqn{z} in Supplement Section 5.1), although the code below allows right_endpt to be a vector
right_endpt_len = length(parameters$right_endpt[1])
group_sub = c(1, 2, 3, 4)
n_agrid_eval = length(as_eval)
# computing the confidence bands
for (j in 1:4) { 
  set.seed(916 + j, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection") # because there's bootstrap
  # find the index for \eqn{\hat{r}} in terms of the \code{as_eval} and \code{as} 
  mu_a0 = apply(Ta[[j]], 2, mean)
  b_a_vec = t(sapply(1:n_agrid_eval, FUN = function (i) {
    min(which(as >= as_eval[i]))
  }))
  right_ind_keep = sapply(right_endpt, function(i) {
    if (sum(mu_a0[b_a_vec] / parameters$tau > i) != 0) {
      max(which(mu_a0[b_a_vec] / parameters$tau > i))
    } else {1}
  }) 
  right_ind_keep_as = sapply(right_endpt, function(i) {
    if (sum(mu_a0 / parameters$tau > i) != 0) {
      max(which(mu_a0 / parameters$tau > i))
    } else {1}
  }) 
  # computing the EL confidence band
  CB_out = elfband(Ta = Ta[[j]], as = as, as_eval = as_eval, n_boot = parameters$nboot, mu_a = NA, alpha = parameters$alpha_vec, as_right_ind_keep0 = right_ind_keep, as_right_ind_keep_as0 = right_ind_keep_as)
  # computing the MFDbs confidence band mentioned in Section 3
  GeoCB_Bs <- fmean(Y = t(Ta[[j]]), as = as, as_eval = as_eval, n_boot = parameters$nboot, alpha = parameters$alpha_vec[1], btypes = c("Bs"), as_right_ind_keep_as0 = right_ind_keep_as, right_endpt = right_endpt)
  for (as_right_ind in 1:right_endpt_len) {
    # computing the Cao confidence band mentioned in Section 3
    CaoCB <- Caoband(Y = Ta[[j]][, 1:right_ind_keep_as[as_right_ind]], as = as[1:right_ind_keep_as[as_right_ind]], checkpts = as_eval[1:right_ind_keep[as_right_ind]], b_M = parameters$nboot, mean_true = mu_a0[1:right_ind_keep_as[as_right_ind]], alpha = parameters$alpha_vec[1]) # Since we do not know \code{mean_true} in real data analysis, here we use the sample mean as \code{mean_true} and do not use the related output \code{cover_or_not} 
    CB_out$CaoCB[[as_right_ind]] = CaoCB$confidence_bounds
    if (n_agrid_eval - right_ind_keep[as_right_ind] > 0) {
      CB_out$CaoCB_Nair[[as_right_ind]] = rbind(CB_out$CaoCB[[as_right_ind]], matrix(-1, ncol = 2, nrow = n_agrid_eval - right_ind_keep[as_right_ind]))
      CB_out$CaoCB_Nair[[as_right_ind]][-(1:right_ind_keep[as_right_ind]), 1] = (CB_out$CaoCB[[as_right_ind]])[right_ind_keep[as_right_ind], 1]
      CaoCB_as_right = (CB_out$CaoCB_Nair[[as_right_ind]])[right_ind_keep[as_right_ind], 2]
      CB_out$CaoCB_Nair[[as_right_ind]][-(1:right_ind_keep[as_right_ind]), 2] = CaoCB_as_right * (CaoCB_as_right < 0)
    } else {
      CB_out$CaoCB_Nair[[as_right_ind]] = CB_out$CaoCB[[as_right_ind]]
    } 
    # computing the Cao2 confidence band mentioned in Section 3
    CaoCB2<- Caoband(Y = Ta[[j]][, 1:right_ind_keep_as[as_right_ind]], as = as[1:right_ind_keep_as[as_right_ind]], checkpts = as_eval[1:right_ind_keep[as_right_ind]], b_M = parameters$nboot, mean_true = mu_a0[1:right_ind_keep_as[as_right_ind]], alpha = parameters$alpha_vec[1], N_m = right_ind_keep_as[as_right_ind] - 2) # Since we do not know \code{mean_true} in real data analysis, here we use the sample mean as \code{mean_true} and do not use the related output \code{cover_or_not} 
    CB_out$CaoCB2[[as_right_ind]]= CaoCB2$confidence_bounds
    if (n_agrid_eval - right_ind_keep[as_right_ind] > 0) {
      CB_out$CaoCB2_Nair[[as_right_ind]] = rbind(CB_out$CaoCB2[[as_right_ind]], matrix(-1, ncol = 2, nrow = n_agrid_eval - right_ind_keep[as_right_ind]))
      CB_out$CaoCB2_Nair[[as_right_ind]][-(1:right_ind_keep[as_right_ind]), 1] = (CB_out$CaoCB2[[as_right_ind]])[right_ind_keep[as_right_ind], 1]
      CaoCB2_as_right = CB_out$CaoCB2_Nair[[as_right_ind]][right_ind_keep[as_right_ind], 2]
      CB_out$CaoCB2_Nair[[as_right_ind]][-(1:right_ind_keep[as_right_ind]), 2] = CaoCB2_as_right * (CaoCB2_as_right < 0)
    } else {
      CB_out$CaoCB2_Nair[[as_right_ind]] = CB_out$CaoCB2[[as_right_ind]]
    } 
    # computing the MFD confidence band mentioned in Section 3
    h <- cv.select(as[1:right_ind_keep_as[as_right_ind]], Ta[[j]][, 1:right_ind_keep_as[as_right_ind]], 1) 
    scbTa <- scb.mean00(as[1:right_ind_keep_as[as_right_ind]], Ta[[j]][, 1:right_ind_keep_as[as_right_ind]], bandwidth = h, degree = 1, gridsize = right_ind_keep[as_right_ind]) 
    CB_out$MFD_CB[[as_right_ind]] = t(scbTa$normscb) 
    if (n_agrid_eval - right_ind_keep[as_right_ind] > 0) {
      CB_out$MFD_CB_Nair[[as_right_ind]] = cbind(CB_out$MFD_CB[[as_right_ind]], matrix(-1, nrow = 2, ncol = n_agrid_eval - right_ind_keep[as_right_ind]))
      CB_out$MFD_CB_Nair[[as_right_ind]][2, -(1:right_ind_keep[as_right_ind])] = (CB_out$MFD_CB[[as_right_ind]])[2, right_ind_keep[as_right_ind]]
      MFD_CB_as_right = CB_out$MFD_CB_Nair[[as_right_ind]][1, right_ind_keep[as_right_ind]]
      CB_out$MFD_CB_Nair[[as_right_ind]][1, -(1:right_ind_keep[as_right_ind])] = MFD_CB_as_right * (MFD_CB_as_right < 0)
    } else {
      CB_out$MFD_CB_Nair[[as_right_ind]] = CB_out$MFD_CB[[as_right_ind]]
    } 
    # computing the Geo confidence band mentioned in Section 3 (\code{btypes = "BEc"} below; \code{btypes = "naive.t"} below is not reported in Section 3 but for comparison purpose only, because it only guarantees pointwise coverage accuracy instead of simultaneous coverage accuracy)
    GeoCB_noBs <- fmean(Y = t(Ta[[j]][, 1:right_ind_keep_as[as_right_ind]]), as = as[1:right_ind_keep_as[as_right_ind]], as_eval = as_eval[1:right_ind_keep[as_right_ind]], n_boot = parameters$nboot, alpha = parameters$alpha_vec[1], btypes = c("BEc", "naive.t"), as_right_ind_keep_as0 = right_ind_keep_as, right_endpt = right_endpt)
    CB_out$GeoCB_noBs[[as_right_ind]] = GeoCB_noBs$confidence_bounds 
    if (n_agrid_eval - right_ind_keep[as_right_ind] > 0) {
      CB_out$GeoCB_Nair_noBs[[as_right_ind]] = rbind(CB_out$GeoCB_noBs[[as_right_ind]], matrix(-1, ncol = 4, nrow = n_agrid_eval - right_ind_keep[as_right_ind]))
      CB_out$GeoCB_Nair_noBs[[as_right_ind]][-(1:right_ind_keep[as_right_ind]), c(1, 3)] = matrix(CB_out$GeoCB_noBs[[as_right_ind]][right_ind_keep[as_right_ind], c(1, 3)], byrow = T, nrow = n_agrid_eval - right_ind_keep[as_right_ind], ncol = 2)
      GeoCB_as_right = CB_out$GeoCB_Nair_noBs[[as_right_ind]][right_ind_keep[as_right_ind], c(2, 4)]
      CB_out$GeoCB_Nair_noBs[[as_right_ind]][-(1:right_ind_keep[as_right_ind]), c(2, 4)] = matrix(GeoCB_as_right * (GeoCB_as_right < 0), byrow = T, nrow = n_agrid_eval - right_ind_keep[as_right_ind], ncol = 2)
    } else {
      CB_out$GeoCB_Nair_noBs[[as_right_ind]] = CB_out$GeoCB_noBs[[as_right_ind]]
    } 
    # continue computing the MFDbs confidence band mentioned in Section 3, based on GeoCB_Bs computed far above (before the for (as_right_ind... loop)
    CB_out$GeoCB_Bs[[as_right_ind]] = GeoCB_Bs$confidence_bounds[, ((2 * as_right_ind - 1):(2 * as_right_ind))]  
    if (n_agrid_eval - right_ind_keep[as_right_ind] > 0) {
      CB_out$GeoCB_Nair_Bs[[as_right_ind]] = CB_out$GeoCB_Bs[[as_right_ind]]
      CB_out$GeoCB_Nair_Bs[[as_right_ind]][-(1:right_ind_keep[as_right_ind]), 1] = rep(CB_out$GeoCB_Bs[[as_right_ind]][right_ind_keep[as_right_ind], 1], times = n_agrid_eval - right_ind_keep[as_right_ind])
      GeoCB_as_right = CB_out$GeoCB_Nair_Bs[[as_right_ind]][right_ind_keep[as_right_ind], 2]
      CB_out$GeoCB_Nair_Bs[[as_right_ind]][-(1:right_ind_keep[as_right_ind]), 2] = rep(GeoCB_as_right * (GeoCB_as_right < 0), times = n_agrid_eval - right_ind_keep[as_right_ind])
    } else {
      CB_out$GeoCB_Nair_Bs[[as_right_ind]] = CB_out$GeoCB_Bs[[as_right_ind]]
    }
  }
  setwd(parameters$dir_path2)
  # saving the resulting file
  save(CB_out,  file = paste("CBgroup_", j, "_group_sub_", paste(group_sub, collapse = "_"), ".Rdata", sep = ""))
}
# reading the resulting file:
CB_out01_mat = array(0, c(4, 9, right_endpt_len)) # number of groups x number of CB types x number of endpoints
CB_lengths_mean_mat = array(0, c(4, 9, right_endpt_len))
CB_mono_mat = array(0, c(4, 9, right_endpt_len))
ELrangej = matrix(0, nrow = 4, ncol = 2)
for (j in 1:4) {
  load(paste("CBgroup_", j, "_group_sub_", paste(group_sub, collapse = "_"), ".Rdata", sep = ""))
  ELrangej[j, ] = range(CB_out$EL_CB_Nair[[1]])
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
    # computing range-violation results of the confidence bands mentioned in Section 3.1
    EL_CB_Nair_out01[as_right_ind] = (sum((((CB_out$EL_CB_Nair[[as_right_ind]])[2, ] > parameters$tau) + ((CB_out$EL_CB_Nair[[as_right_ind]])[1, ] < 0)) > 0) >= 1)  
    EP_CB_Nair_out01[as_right_ind] = (sum((((CB_out$EP_CB_Nair[[as_right_ind]])[2, ] > parameters$tau) + ((CB_out$EP_CB_Nair[[as_right_ind]])[1, ] < 0)) > 0) >= 1)  
    HW_CB_out01[as_right_ind] = (sum((((CB_out$HW_CB_eval)[2, , as_right_ind] > parameters$tau) + ((CB_out$HW_CB_eval)[1, , as_right_ind] < 0)) > 0) >= 1)  
    Cao_CB_out01[as_right_ind] = (sum((((CB_out$CaoCB_Nair[[as_right_ind]])[, 1] > parameters$tau) + ((CB_out$CaoCB_Nair[[as_right_ind]])[, 2] < 0)) > 0) >= 1)  
    Cao_CB2_out01[as_right_ind] = (sum((((CB_out$CaoCB2_Nair[[as_right_ind]])[, 1] > parameters$tau) + ((CB_out$CaoCB2_Nair[[as_right_ind]])[, 2] < 0)) > 0) >= 1)  
    Bs_CB_out01[as_right_ind] = (sum((((CB_out$GeoCB_Nair_Bs[[as_right_ind]])[, 1] > parameters$tau) + ((CB_out$GeoCB_Nair_Bs[[as_right_ind]])[, 2] < 0)) > 0) >= 1)  
    BEc_CB_out01[as_right_ind] = (sum((((CB_out$GeoCB_Nair_noBs[[as_right_ind]])[, 1] > parameters$tau) + ((CB_out$GeoCB_Nair_noBs[[as_right_ind]])[, 2] < 0)) > 0) >= 1)  
    naivet_CB_out01[as_right_ind] = (sum((((CB_out$GeoCB_Nair_noBs[[as_right_ind]])[, 3] > parameters$tau) + ((CB_out$GeoCB_Nair_noBs[[as_right_ind]])[, 4] < 0)) > 0) >= 1)  
    MFD_CB_out01[as_right_ind] = (sum((((CB_out$MFD_CB_Nair[[as_right_ind]])[2, ] > parameters$tau) + ((CB_out$MFD_CB_Nair[[as_right_ind]])[1, ] < 0)) > 0) >= 1)  
    CB_out01_mat[j, , as_right_ind] = c(EL_CB_Nair_out01[as_right_ind], EP_CB_Nair_out01[as_right_ind], HW_CB_out01[as_right_ind], Cao_CB_out01[as_right_ind], Cao_CB2_out01[as_right_ind], Bs_CB_out01[as_right_ind], BEc_CB_out01[as_right_ind], naivet_CB_out01[as_right_ind], MFD_CB_out01[as_right_ind]) 
    # computing average widths of the confidence bands mentioned in Section 3.1
    EL_CB_Nair_length[as_right_ind] = mean((CB_out$EL_CB_Nair[[as_right_ind]])[2, ] - (CB_out$EL_CB_Nair[[as_right_ind]])[1, ])
    EP_CB_Nair_length[as_right_ind] = mean((CB_out$EP_CB_Nair[[as_right_ind]])[2, ] - (CB_out$EP_CB_Nair[[as_right_ind]])[1, ])
    HW_CB_length[as_right_ind] = mean((CB_out$HW_CB_eval)[2, , as_right_ind] - (CB_out$HW_CB_eval)[1, , as_right_ind])
    Cao_CB_length[as_right_ind] = mean((CB_out$CaoCB_Nair[[as_right_ind]])[, 1] - (CB_out$CaoCB_Nair[[as_right_ind]])[, 2])
    Cao_CB2_length[as_right_ind] = mean((CB_out$CaoCB2_Nair[[as_right_ind]])[, 1] - (CB_out$CaoCB2_Nair[[as_right_ind]])[, 2])
    Bs_CB_length[as_right_ind] = mean((CB_out$GeoCB_Nair_Bs[[as_right_ind]])[, 1] - (CB_out$GeoCB_Nair_Bs[[as_right_ind]])[, 2])
    BEc_CB_length[as_right_ind] = mean((CB_out$GeoCB_Nair_noBs[[as_right_ind]])[, 1] - (CB_out$GeoCB_Nair_noBs[[as_right_ind]])[, 2])
    naivet_CB_length[as_right_ind] = mean((CB_out$GeoCB_Nair_noBs[[as_right_ind]])[, 3] - (CB_out$GeoCB_Nair_noBs[[as_right_ind]])[, 4])
    MFD_CB_length[as_right_ind] = mean((CB_out$MFD_CB_Nair[[as_right_ind]])[2, ] - (CB_out$MFD_CB_Nair[[as_right_ind]])[1, ])
    CB_lengths_mean_mat[j, , as_right_ind] = c(EL_CB_Nair_length[as_right_ind], EP_CB_Nair_length[as_right_ind], HW_CB_length[as_right_ind], Cao_CB_length[as_right_ind], Cao_CB2_length[as_right_ind], Bs_CB_length[as_right_ind], BEc_CB_length[as_right_ind], naivet_CB_length[as_right_ind], MFD_CB_length[as_right_ind]) 
    # computing monotonicity preservation of the confidence bands mentioned in Section 3.1
    EL_Nair_mono[as_right_ind] = (sum(diff((CB_out$EL_CB_Nair[[as_right_ind]])[2, ] ) > 0) == 0) + (sum(diff((CB_out$EL_CB_Nair[[as_right_ind]])[1, ] ) > 0) == 0)
    EP_Nair_mono[as_right_ind] = (sum(diff((CB_out$EP_CB_Nair[[as_right_ind]])[2, ] ) > 0) == 0) + (sum(diff((CB_out$EP_CB_Nair[[as_right_ind]])[1, ] ) > 0) == 0)
    HW_mono[as_right_ind] = (sum(diff((CB_out$HW_CB_eval)[2, , as_right_ind]) > 0) == 0) + (sum(diff((CB_out$HW_CB_eval)[1, , as_right_ind]) > 0) == 0)
    Cao_mono[as_right_ind] = (sum(diff((CB_out$Cao_CB_Nair[[as_right_ind]])[, 1]) > 0) == 0) + (sum(diff((CB_out$Cao_CB_Nair[[as_right_ind]])[, 2]) > 0) == 0)
    Cao2_mono[as_right_ind] = (sum(diff((CB_out$Cao_CB2_Nair[[as_right_ind]])[, 1]) > 0) == 0) + (sum(diff((CB_out$Cao_CB2_Nair[[as_right_ind]])[, 2]) > 0) == 0)
    Bs_mono[as_right_ind] = (sum(diff((CB_out$GeoCB_Nair_Bs[[as_right_ind]])[, 1]) > 0) == 0) + (sum(diff((CB_out$GeoCB_Nair_Bs[[as_right_ind]])[, 2]) > 0) == 0)
    BEc_mono[as_right_ind] = (sum(diff((CB_out$GeoCB_Nair_noBs[[as_right_ind]])[, 1]) > 0) == 0) + (sum(diff((CB_out$GeoCB_Nair_noBs[[as_right_ind]])[, 2]) > 0) == 0)
    naivet_mono[as_right_ind] = (sum(diff((CB_out$GeoCB_Nair_noBs[[as_right_ind]])[, 3]) > 0) == 0) + (sum(diff((CB_out$GeoCB_Nair_noBs[[as_right_ind]])[, 4]) > 0) == 0)
    MFD_mono[as_right_ind] = (sum(diff((CB_out$MFD_CB_Nair[[as_right_ind]])[2, ]) > 0) == 0) + (sum(diff((CB_out$MFD_CB_Nair[[as_right_ind]])[1, ]) > 0) == 0)
    CB_mono_mat[j, , as_right_ind] = c(EL_Nair_mono[as_right_ind], EP_Nair_mono[as_right_ind], HW_mono[as_right_ind], Cao_mono[as_right_ind], Cao2_mono[as_right_ind], Bs_mono[as_right_ind], BEc_mono[as_right_ind], naivet_mono[as_right_ind], MFD_mono[as_right_ind]) 
  }
}
dimnames(CB_out01_mat)[[2]] <- c('EL', 'EP', 'NS', 'Cao', 'Cao2', 'MFDbs', 'Geo', 'naivet', 'MFD')
dimnames(CB_lengths_mean_mat)[[2]] <- c('EL', 'EP', 'NS', 'Cao', 'Cao2', 'MFDbs', 'Geo', 'naivet', 'MFD')
dimnames(CB_mono_mat)[[2]] <- c('EL', 'EP', 'NS', 'Cao', 'Cao2', 'MFDbs', 'Geo', 'naivet', 'MFD')
# ------------ Figure 5 -------------- #
setwd(parameters$dir_path2)
left_mar=6
bot_mar=4
n_quant = 5
time_unit = 24
text_line = 2.5
lwds = 2.2
png("paper_fnl_fig_NHANES_CBs_20220104_test.png",width = 16, height = 6, units='in', res = 600)
require(graphics)
par(mfrow = c(1,3),mar = c(bot_mar, left_mar, 1, 1.1), oma = c(0.3, 0.5, 0, 0), las = 1, mgp = c(3, 1, 0), cex.axis = 1.5, font.lab = 50, tck = -0.015) 
# left panel of Figure 5
plot(as, (Ta[[1]])[1, ], type = "n", col = 'gray', ylim = c(min(unlist(Ta)), max(unlist(Ta))) * time_unit, xlab = '', ylab = '')
mtext(expression(T[j](a)), cex = 1.3, line = text_line, side = 2)
mtext('activity level (a)', cex = 1.3, line = 3, side = 1)
j = 3
rep_subject = 1:(n_quant + 1) * 0
for (quanti in 0:n_quant) {
  quant = quantile((Ta[[j]])[, 1], probs = quanti * (1 / n_quant))
  rep_subject[quanti + 1] = which.min(abs((Ta[[j]])[, 1] - quant))
  points(as, (Ta[[j]])[rep_subject[quanti + 1], ] * time_unit, type = "l", col = "red", lwd = lwds - 0.3) 
}
j = 2
rep_subject = 1:(n_quant + 1) * 0
for (quanti in 0:n_quant) {
  quant = quantile((Ta[[j]])[, 1], probs=quanti * (1 / n_quant))
  rep_subject[quanti + 1] = which.min(abs((Ta[[j]])[, 1] - quant))
  points(as, (Ta[[j]])[rep_subject[quanti + 1], ] * time_unit, type = "l", col = "green", lwd = lwds - 0.3)  
}
j = 1
rep_subject = 1:(n_quant + 1) * 0
for (quanti in 0:n_quant) {
  quant = quantile((Ta[[j]])[, 1], probs = quanti * (1 / n_quant))
  rep_subject[quanti + 1] = which.min(abs((Ta[[j]])[, 1] - quant))
  points(as, (Ta[[j]])[rep_subject[quanti + 1], ] * time_unit, type = "l", lwd = lwds - 0.3) 
}
j = 4
rep_subject = 1:(n_quant + 1) * 0
for (quanti in 0:n_quant) {
  quant = quantile((Ta[[j]])[, 1], probs = quanti * (1 / n_quant))
  rep_subject[quanti + 1] = which.min(abs((Ta[[j]])[, 1] - quant))
  points(as, (Ta[[j]])[rep_subject[quanti + 1], ] * time_unit, type = "l", col = "purple", lwd = lwds + 1, lty = 3) 
}
# middle panel of Figure 5
plot(as, (CB_out$EL_CB_Nair[[as_right_ind]])[1, ] * time_unit, type = "n", ylim = c(min(ELrangej), max(ELrangej)) * time_unit, xlab = '', ylab = '')
mtext(expression(mu[j](a)),cex = 1.3, line = text_line, side=2)
mtext('activity level (a)',cex = 1.3, line = 3, side=1)
as_right_ind = 1
j = 1
load(paste("CBgroup_", j, "_group_sub_", paste(group_sub, collapse = "_"), ".Rdata", sep = ""))
points(as, CB_out$mu_hat_vec * time_unit, type="l", lty = 2, lwd = lwds)
points(as, (CB_out$EL_CB_Nair[[as_right_ind]])[1, ] * time_unit, type="l", lwd = lwds)
points(as, (CB_out$EL_CB_Nair[[as_right_ind]])[2, ] * time_unit, type="l", lwd = lwds)
j = 3
load(paste("CBgroup_", j, "_group_sub_", paste(group_sub, collapse = "_"), ".Rdata", sep = ""))
points(as, CB_out$mu_hat_vec * time_unit, type = "l", col = "red", lty = 2, lwd = lwds)
points(as, (CB_out$EL_CB_Nair[[as_right_ind]])[1, ] * time_unit, type = "l", col = "red", lwd = lwds) 
points(as, (CB_out$EL_CB_Nair[[as_right_ind]])[2, ] * time_unit, type = "l", col = "red", lwd = lwds)
# right panel of Figure 5
xlim_indx = 91:101
xlims = c(90.368, 99.63)
gaps = c(3.74, 4.205)  
j = 1
load(paste("CBgroup_", j, "_group_sub_", paste(group_sub, collapse = "_"), ".Rdata", sep = ""))
par(bty="o") 
gap.plot(as[xlim_indx ], (CB_out$EL_CB_Nair[[as_right_ind]])[1, xlim_indx ] * time_unit, gap = gaps, gap.axis = "y", ytics = seq(3.6, 4.4, by = 0.1), yticlab = seq(3.6, 4.4, by = 0.1), type = "n", ylim = c(0.15, 0.185) * time_unit, xlim = xlims, xlab = '', ylab = '')
abline(h = seq(3.74, 3.749, by = 0.0001), col = "white")  # hiding horizontal lines from \code{gap.plot}
axis.break(axis = 2, breakpos = 3.745)
abline(v=100.004, lwd = 1.5)  # completing missing vertical line due to \code{gap.plot}
mtext('activity level (a)', cex = 1.3, line = 3, side = 1)
mtext(expression(mu[1](a)), cex = 1.3, line = text_line, side = 2, las = 2)
# plotting the EL confidence band mentioned in Section 3
gap.plot(as[xlim_indx], (CB_out$EL_CB_Nair[[as_right_ind]])[1, xlim_indx] * time_unit, lwd = lwds, type = "l", lty = 2, add = TRUE, gap = gaps, gap.axis = "y", ylim = c(0.15, 0.185) * time_unit, xlim = xlims)
gap.plot(as[xlim_indx], (CB_out$EL_CB_Nair[[as_right_ind]])[2, xlim_indx] * time_unit, lwd = lwds, type = "l", lty = 2, add = TRUE, gap = gaps, gap.axis = "y", ylim = c(0.15, 0.185) * time_unit, xlim = xlims)
# plotting the EP confidence band mentioned in Section 3
gap.plot(as[xlim_indx], (CB_out$EP_CB_Nair[[as_right_ind]])[1, xlim_indx] * time_unit, type = "l", lwd = lwds, col = "#FF99FF", add = TRUE, gap = gaps, gap.axis="y", ylim = c(0.15, 0.185) * time_unit, xlim = xlims)  
gap.plot(as[xlim_indx], (CB_out$EP_CB_Nair[[as_right_ind]])[2, xlim_indx] * time_unit, type = "l", lwd = lwds, col = "#FF99FF", add = TRUE, gap = gaps, gap.axis="y", ylim = c(0.15, 0.185) * time_unit, xlim = xlims)
# plotting the NS confidence band mentioned in Section 3
gap.plot(as[xlim_indx], (CB_out$HW_CB_eval)[1, xlim_indx, as_right_ind] * time_unit, type = "l", col = "purple", add = TRUE, gap = gaps, gap.axis = "y", ylim = c(0.15, 0.185) * time_unit, xlim = xlims) 
gap.plot(as[xlim_indx], (CB_out$HW_CB_eval)[2, xlim_indx, as_right_ind] * time_unit, type = "l", col = "purple", add = TRUE, gap = gaps, gap.axis = "y", ylim = c(0.15, 0.185) * time_unit, xlim = xlims)
# plotting the MFDbs confidence band mentioned in Section 3
gap.plot(as[xlim_indx], (CB_out$GeoCB_Nair_Bs[[as_right_ind]])[xlim_indx, 2] * time_unit, type = "l", lty = 1, lwd = lwds, col = "orange", add = TRUE, gap = gaps, gap.axis = "y", ylim = c(0.15, 0.185) * time_unit, xlim = xlims)
gap.plot(as[xlim_indx], (CB_out$GeoCB_Nair_Bs[[as_right_ind]])[xlim_indx, 1] * time_unit, type = "l", lty = 1, lwd = lwds, col = "orange", add = TRUE, gap = gaps, gap.axis = "y", ylim = c(0.15, 0.185) * time_unit, xlim = xlims)
# plotting the Geo confidence band mentioned in Section 3
gap.plot(as[xlim_indx], (CB_out$GeoCB_Nair_noBs[[as_right_ind]])[xlim_indx, 2] * time_unit, type = "l", lty = 1, lwd = lwds, col = "lightgreen", add = TRUE, gap = gaps, gap.axis = "y", ylim = c(0.15, 0.185) * time_unit, xlim = xlims)
gap.plot(as[xlim_indx], (CB_out$GeoCB_Nair_noBs[[as_right_ind]])[xlim_indx, 1] * time_unit, type = "l", lty = 1, lwd = lwds, col = "lightgreen", add = TRUE, gap = gaps, gap.axis = "y", ylim = c(0.15, 0.185) * time_unit, xlim = xlims)
dev.off()
# ------------ remaining results in Section 4 -------------- #
# result 6 in Section 4---comparing the average widths of the EL, EP and NS confidence bands
CB_lengths_mean_mat[c(1, 3), 1:3,1]
# result 7 in Section 4---comparing the widths of the EL, MFDbs and Geo confidence bands over the entire ragne of activity levels of interest
j=1
as_right_ind = 1
load(paste("CBgroup_", j, "_group_sub_", paste(group_sub, collapse = "_"), ".Rdata", sep = ""))
EL_CB_Nair_widths = (CB_out$EL_CB_Nair[[as_right_ind]])[2, ] - (CB_out$EL_CB_Nair[[as_right_ind]])[1, ]
Bs_CB_widths = (CB_out$GeoCB_Nair_Bs[[as_right_ind]])[, 1] - (CB_out$GeoCB_Nair_Bs[[as_right_ind]])[, 2]
BEc_CB_widths = (CB_out$GeoCB_Nair_noBs[[as_right_ind]])[, 1] - (CB_out$GeoCB_Nair_noBs[[as_right_ind]])[, 2]
length(EL_CB_Nair_widths) 
sum(Bs_CB_widths < EL_CB_Nair_widths) # the MFDbs band is wider than EL at all points
which(BEc_CB_widths  < EL_CB_Nair_widths)
sort(BEc_CB_widths - EL_CB_Nair_widths, index.return = TRUE, decreasing = TRUE)$ix # the Geo band is wider than EL in the extremes of the activity levels
# result 8 in Section 4---all the bands not based on smoothing respect range and monotonicity constraints of the occupation time data
CB_out01_mat[c(1, 3), -c(4:5, 8:9), 1] # for subgroups 1 (first row) and 3 (second row), the range-violation results of the confidence bands that are not based on smoothing mentioned in Section 4
CB_mono_mat[c(1, 3), -c(4:5, 8:9), 1] # for subgroups 1 (first row) and 3 (second row), the monotonicity preservation results of the confidence bands that are not based on smoothing mentioned in Section 4; see the caption of Table 1 for what the number means
# result 9 in Section 4---the 95% EL confidence interval for the mean hours per day of sedentary behavior, for cut-points 50, 100, 500, respectively
(a_indx = c(which.min(abs(as-49)), which.min(abs(as-99)), which.min(abs(as-499)))) # Since we are interested in time spent < 50, 100 or 500 counts/min, or equivalently <= 49, 99 or 499 counts/min, here we find the indices corresponding to 49, 99 and 499 counts/min, respectively.
j = 1
as_right_ind = 1
load(paste("CBgroup_", j, "_group_sub_", paste(group_sub, collapse = "_"), ".Rdata", sep = ""))
(1 - (CB_out$EL_CB_Nair[[as_right_ind]])[2, a_indx]) * 24 # the lower bounds of the EL confidence interval for the mean hours per day of sedentary behavior, for cut-points 50, 100, 500, respectively
(1 - (CB_out$EL_CB_Nair[[as_right_ind]])[1, a_indx]) * 24 # the upper bounds of the EL confidence interval for the mean hours per day of sedentary behavior, for cut-points 50, 100, 500, respectively
