product_mat <- function(x, y) {
  # A modified * function so that  0 * Inf = 0
  # Args:
  #   x, y: two numbers or vectors to be multiplied
  # Returns:
  #   x * y, where an element becomes 0 if the corresponding element of one of x or y is 0
  out <- x*y
  for (i in 1:dim(x)[1]) {
      out[i, (x[i, ] == 0 | y[i, ] == 0)] <- 0
  }
  return (out)
}
get.crit.supnorm.simple = function (cov.m, n.sim, prob, as_right_ind_keep_as0) {
  # Given the estimated functional mean and variance, compute the Geo confidence bands over \code{as} (defined in the function \code{fmean}) modified from the bands by Choi and Reimherr (2018) mentioned in Section 3.
  # Args:
  #   cov.m: a vector containing the estimated functional variance over \code{as} (defined in the function \code{fmean})
  #   n.sim: number of bootstrap samples involved in the \code{"Bs"} band.
  #   prob: the given probability that the critical value corresponds to; typically \code{1-alpha} (\code{alpha} is defined in the function \code{fmean})
  # Returns:
  #   a vector containing the critical values, with the (\code{i})-th element corresponding to \code{right_endpt[i]} (defined in the function \code{fregion.band})
  # Authors:
  #   The code is directly modified by Hsin-wen Chang based on the \code{get.crit.supnorm.simple} function in the the R package \pkg{fregion} (version 0.0934) by Hyunphil Choi
  # References:
  #   Choi, H. and Reimherr, M. (2018) A geometric approach to confidence regions and bands for functional parameters. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 80, 239--260.
  # Source:
  #   The R package \pkg{fregion} (version 0.0934) from https://github.com/hpchoi/fregion
  cov.m_full = cov.m
  J <- length(diag(cov.m))
  simu = matrix(rnorm(J * n.sim), ncol = n.sim)
  right_endpt_len = length(as_right_ind_keep_as0)
  out = 1:right_endpt_len * 0
  for (as_right_ind in 1:right_endpt_len) {
    cov.m = cov.m_full[1:as_right_ind_keep_as0[as_right_ind], 1:as_right_ind_keep_as0[as_right_ind]]
    J <- length(diag(cov.m))
    X <- array(0, dim = c(J, n.sim))
    sd.v <- sqrt(diag(cov.m))
    cor.m <- product_mat(cov.m, 1 / outer(sd.v, sd.v))
    eigen <- eigen(cor.m)
    eigen$values[eigen$values < .Machine$double.eps] <- 0
    sd.m.2 <- crossprod(t(eigen$vectors), diag(sqrt(eigen$values)))
    X <- crossprod(t(sd.m.2), simu[1:as_right_ind_keep_as0[as_right_ind], ])
    MaxX <- apply(abs(X), 2, max)
    out[as_right_ind] = quantile(MaxX, prob)
  }
  return(out)
}
make.band.Bs = function (cov, conf.level, sim.size = n_boot, as_right_ind_keep_as0) {
  # Given the estimated functional mean and variance, compute the Geo confidence bands over \code{as} (defined in the function \code{fmean}) modified from the bands by Choi and Reimherr (2018) mentioned in Section 3.
  # Args:
  #   cov: a vector containing the estimated functional variance over \code{as} (defined in the function \code{fmean})
  #   conf.level: \code{1-alpha} (\code{alpha} is defined in the function \code{fmean})
  #   sim.size: number of bootstrap samples involved in the \code{"Bs"} band.
  #   as_right_ind_keep_as0: This is the index in the vector \code{as} (defined in the function \code{fmean}) corresponding to \eqn{\hat{r}} defined in Supplement Section 5.1, and it can be a vector that corresponds to different values of \eqn{z} in Supplement Section 5.1.
  # Returns:
  #   a matrix containing half the width of the \code{"Bs"} band, with the (\code{i})-th column corresponding to \code{right_endpt[i]} (defined in the function \code{fregion.band})
  # Authors:
  #   The code is directly modified by Hsin-wen Chang based on the \code{make.band.Bs} function in the the R package \pkg{fregion} (version 0.0934) by Hyunphil Choi
  # References:
  #   Choi, H. and Reimherr, M. (2018) A geometric approach to confidence regions and bands for functional parameters. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 80, 239--260.
  # Source:
  #   The R package \pkg{fregion} (version 0.0934) from https://github.com/hpchoi/fregion
  cov.m <- cov
  crit.Bs <- get.crit.supnorm.simple(cov.m = cov.m, n.sim = sim.size, prob = conf.level, as_right_ind_keep_as0 = as_right_ind_keep_as0)
  band.eval <- sqrt(diag(cov.m)) %o% crit.Bs
  return(band.eval)
}
fregion.band <- function(x, cov, N = 1, type = c("Bs", "BEc"), conf.level = c(0.95), sim.size = 10000, as_right_ind_keep_as0, right_endpt) {
  # Given the estimated functional mean and variance, compute the Geo confidence bands over \code{as} (defined in the function \code{fmean}) modified from the bands by Choi and Reimherr (2018) mentioned in Section 3.
  # Args:
  #   x: a vector containing the estimated functional mean over \code{as} (defined in the function \code{fmean})
  #   cov: a vector containing the estimated functional variance over \code{as} (defined in the function \code{fmean})
  #   N: the number of observed functional curves or subjects
  #   type: the type(s) of bands used in fregion.band: \code{"BEc"}, \code{"Bs"}, or \code{"naive.t"}. (default = \code{"Bs"}) The outputs from the \code{"BEc"} and \code{"naive.t"} bands are relatively simple, so they can be requested together via \code{btypes = c("BEc", "naive.t")}. On the contrary, the \code{"Bs"} band can only be requested alone, because the \code{"Bs"} band utilizes the Nair's two-step approach described in Supplement Section 5.1, which makes the length of the output of the \code{"Bs"} band depends on the length of \code{right_endpt}.
  #   conf.level: \code{1-alpha} (\code{alpha} is defined in the function \code{fmean})
  #   sim.size: number of bootstrap samples involved in the \code{"Bs"} band.
  #   as_right_ind_keep_as0: This is the index in the vector \code{as} (defined in the function \code{fmean}) corresponding to \eqn{\hat{r}} defined in Supplement Section 5.1, and it can be a vector that corresponds to different values of \eqn{z} in Supplement Section 5.1.
  #   right_endpt: the value(s) of \eqn{z} in Supplement Section 5.1.
  # Returns:
  #   a matrix containing the estimated functional mean in the first column and the confidence bounds over \code{as} (defined in the function \code{fmean}) in the later columns. If \code{type == "Bs"}, the (\code{2 * i})-th and (\code{2 * i + 1})-th columns contain the upper and lower bound of the \code{"Bs"} confidence band corresponding to \code{right_endpt[i]} (defined above), respectively. If \code{type != "Bs"}, the (\code{2 * i})-th and (\code{2 * i + 1})-th columns contain the upper and lower bound of the \code{type[i]} (defined above) confidence band, respectively.
  # Authors:
  #   The code is directly modified by Hsin-wen Chang based on the \code{fregion.band} function in the the R package \pkg{fregion} (version 0.0934) by Hyunphil Choi
  # References:
  #   Choi, H. and Reimherr, M. (2018) A geometric approach to confidence regions and bands for functional parameters. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 80, 239--260.
  # Source:
  #   The R package \pkg{fregion} (version 0.0934) from https://github.com/hpchoi/fregion
  ### 1. Check the data type ###
  if ((inherits(x, "numeric") | inherits(x, "matrix"))  & (inherits(cov, "matrix") | inherits(cov, "list") | inherits(cov, "eigen") )) {
    datatype = "vector"
  } else stop ("The format of data is unknown")
  ### 2. Evaluate x and cov ###
  # evaluate x and cov if datatype is "fd".
  # Since all functions for generating bands evaluate fd object inside, we do this here and just use vector/matrix version
  cov.m <- cov
  x.v <- x
  p <- dim(cov.m)[1]
  if (!isSymmetric(cov.m)) cov.m <- (cov.m + t(cov.m)) / 2  # force cov.m to be symmetric
  ## 3. Take loop for conf.level
  result <- as.matrix(x.v, ncol = 1) ;
  colnames(result) <- c("x")
  ## 4. Take eigen decomposition if BEc or BEc1 is used.
  if (  (sum(c("BEc", "BEc1") %in% type) > 0) | (length(grep("BEPC.", type)) > 0) ) {
    eigen.cov.m <- eigen(cov.m) ; eigen.cov.m$values[ eigen.cov.m$values < 0 ] <- 0
  }
  for (i in c(1:length(conf.level))){
    level <- conf.level[i]
    # 5. Find number of fpc to use.
    if (!(level > 0 & level < 1)) stop("conf.level should have values between 0 and 1")
    # 6. Make bands
    if ("Bs" %in% type) {
      right_endpt_len = length(as_right_ind_keep_as0)
      tmp.colnames <- c(colnames(result), c(sapply(1:right_endpt_len, FUN = function (i) {
        c(paste0("Bs.u.", level, "endpt", right_endpt[i]), paste0("Bs.l.", level, "endpt", right_endpt[i]))
      })))
      Bs <- make.band.Bs(cov = cov.m, conf.level = level, sim.size = sim.size, as_right_ind_keep_as0 = as_right_ind_keep_as0) / sqrt(N)
      Bs_bands = do.call(cbind,lapply(1:right_endpt_len, FUN = function (i) {
        cbind(x.v + Bs[,i], x.v - Bs[,i])
      }))
      result <- cbind(result, Bs_bands)
      colnames(result) <- tmp.colnames
    }
    if ("BEc" %in% type) {
      tmp.colnames <- c(colnames(result), paste0("BEc.u.", level), paste0("BEc.l.", level))
      BEc <- fregion::make.band.BEc(eigen = eigen.cov.m, conf.level = level) / sqrt(N)
      result <- cbind(result, x.v + BEc, x.v - BEc);
      colnames(result) <- tmp.colnames
    }
    if ("naive.t" %in% type) {
      tmp.colnames <- c(colnames(result), paste0("naive.t.u.", level), paste0("naive.t.l.", level))
      naive.t <- fregion::make.band.naive.t(cov.m, conf.level = level, df = N - 1) / sqrt(N)
      result <- cbind(result, x.v + naive.t, x.v - naive.t);
      colnames(result) <- tmp.colnames
    }
    tmpJs.location <- grep("BEPC.", type)
    if (length(tmpJs.location) > 0) {
      tmpJs <- as.integer(sub("BEPC.", "", type[tmpJs.location]))
      for (j in tmpJs){
        tmp.colnames <- c(colnames(result), paste0("BEPC.", j, ".u.", level), paste0("BEPC.", j, ".l.", level))
        BEPC <- fregion::make.band.BEPC(eigen = eigen.cov.m, conf.level = level, J = j) / sqrt(N)
        result <- cbind(result, x.v + BEPC, x.v - BEPC);
        colnames(result) <- tmp.colnames
      }
    }
  }
  class(result) <- "fregion.band"
  return(result)
}
fmean <- function(Y, as, as_eval, n_boot = 1000, alpha = 0.05, btypes = "Bs", as_right_ind_keep_as0, right_endpt) {
  # Compute the Geo confidence bands modified from the bands by Choi and Reimherr (2018) mentioned in Section 3
  # Args:
  #   Y: a matrix of data, where each column is a discretized version of a function and each row corresponds to each design time point.
  #   as: a vector consisting of the design time points of the observed functional data; that is, \eqn{{\bf G}_n} defined in Section 2.2 (from the smallest to the largest).
  #   as_eval: a vector consisting of the time points (from the smallest to the largest; e.g., the activity levels of interest) at which we want to evaluate the confidence band. This vector can be denser or coarser than \code{as}, as long as the largest point of \code{as_eval} is no greater than the last point of \code{as} (defined above).
  #   n_boot: number of bootstrap samples for \code{btypes = "Bs"} (see below)
  #   alpha: a number, indicating the significance level of interest (default = 0.05)
  #   btypes: the type(s) of bands used in fregion.band: \code{"BEc"}, \code{"Bs"}, or \code{"naive.t"}. (default = \code{"Bs"}) The outputs from the \code{"BEc"} and \code{"naive.t"} bands are relatively simple, so they can be requested together via \code{btypes = c("BEc", "naive.t")}. On the contrary, the \code{"Bs"} band can only be requested alone, because the \code{"Bs"} band utilizes the Nair's two-step approach described in Supplement Section 5.1, which makes the length of the output of the \code{"Bs"} band depends on the length of \code{right_endpt}.
  #   as_right_ind_keep_as0: This is the index in the vector \code{as} (defined above) corresponding to \eqn{\hat{r}} defined in Supplement Section 5.1, and it can be a vector that corresponds to different values of \eqn{z} in Supplement Section 5.1.
  #   right_endpt: the value(s) of \eqn{z} in Supplement Section 5.1.
  # Returns:
  #   confidence_bounds: a matrix containing the confidence bounds. If \code{btypes == "Bs"}, the (\code{2 * i - 1})-th and (\code{2 * i})-th columns contain the upper and lower bound of the \code{"Bs"} confidence band corresponding to \code{right_endpt[i]}, respectively. If \code{btypes != "Bs"}, the (\code{2 * i - 1})-th and (\code{2 * i})-th columns contain the upper and lower bound of the \code{btypes[i]} confidence band, respectively.
  #   confidence_level: \code{1-alpha}
  #   runtime: the runtime of the \code{fmean} funciton
  # Authors:
  #   The code is first written by Yan-Yu Chen on 2021/03/11 and later modified by Hsin-wen Chang.
  # References:
  #   Choi, H. and Reimherr, M. (2018) A geometric approach to confidence regions and bands for functional parameters. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 80, 239--260.
  # Source:
  #   The R package fregion (version 0.0934) from https://github.com/hpchoi/fregion
  start_time <- Sys.time()
  n_agrid_eval = length(as_eval)
  N <- dim(Y)[2]
  Y_bar <- apply(Y, 1, mean)
  Y_cov <- (Y - Y_bar) %*% t(Y - Y_bar) / (N - 1)
  right_endpt_len = length(as_right_ind_keep_as0)
  Y_bands <-  fregion.band(x = Y_bar, cov = Y_cov, N = N, type = btypes, conf.level = 1-alpha, sim.size = n_boot, as_right_ind_keep_as0 = as_right_ind_keep_as0, right_endpt = right_endpt)
  b_a_vec = t(sapply(1:n_agrid_eval, FUN = function (i) {
    min(which(as >= as_eval[i]))
  }))
  if ("Bs" %in% btypes) {
    nb <- right_endpt_len
    Y_bands.ub <- Y_bands[b_a_vec, c(1:nb) * 2, drop = F]
    Y_bands.lb <- Y_bands[b_a_vec, c(1:nb) * 2 + 1, drop = F]
  } else {
    nb <- length(btypes)
    Y_bands.ub <- Y_bands[b_a_vec, c(1:nb) * 2, drop = F]
    Y_bands.lb <- Y_bands[b_a_vec, c(1:nb) * 2 + 1, drop = F]
  }
  end_time <- Sys.time()
  res <- list('confidence_bounds' = Y_bands[b_a_vec, -1], 'confidence_level' = 1 - alpha, 'runtime' = end_time - start_time)
  return(res)
}
