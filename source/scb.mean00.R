division00 <- function(x, y) {  
  # A modified / function so that  0 / 0 = 0
  # Args:
  #   x, y: two numbers or vectors
  # Returns:
  #   x / y, where an element becomes 0 if the corresponding element of one of x or y is 0
  out <- x/y
  out[as.logical((x == 0)*(x == y))] <- 0
  return (out)
}
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
scb.mean00 = function (x, y, bandwidth, level = 0.95, degree = 1, gridsize = length(x),
                     keep.y = TRUE, nrep = 20000) {
  # Modifying the confidence bands computed by \code{\link[SCBmeanfd:scb.mean]{scb.mean}} (version 1.2.2) to take care of the 0 / 0 case (by setting it to 0). 
  # All the arguments and return values of this modified function are the same as the ones related to the `normal' scbtype in \code{\link[SCBmeanfd:scb.mean]{scb.mean}} (version 1.2.2).
  # Source:
  #   See \href{https://CRAN.R-project.org/package=SCBmeanfd}{R package SCBmeanfd (version 1.2.2)} for the arguments (x, y, bandwidth, level, degree, gridsize, keep.y, nrep) and return values (without qboot, nboot, bootscb).
  caLL <- match.call()
  stopifnot(is.matrix(y))
  stopifnot(all(!is.na(x) & !is.na(y)))
  stopifnot(length(x) == ncol(y))
  n <- nrow(y)
  N <- ncol(y)
  y.hat <- apply(y, 1, function(z) KernSmooth::locpoly(x, z, degree = degree,
                                           bandwidth = bandwidth, gridsize = gridsize)$y)
  mu.hat <- rowMeans(y.hat)
  sigma.hat <- apply(y.hat, 1, sd)
  se <- sigma.hat / sqrt(n)
  r <- y.hat - mu.hat
  lb.norm = ub.norm = q.norm = NULL # modified
  sigma.hat_mat = matrix(rep(sigma.hat, times = dim(r)[2]), nrow = dim(r)[1], ncol = dim(r)[2]) # added
  svd.r <- svd(product_mat(r, 1 / sigma.hat_mat) / sqrt(n - 1), nv = 0) # modified
  ncomp <- which.max(cumsum(svd.r$d ^ 2) > 0.99 * sum(svd.r$d ^ 2))
  vars <- matrix(rnorm(ncomp * nrep), ncomp, nrep)
  M <- svd.r$u[, 1:ncomp] %*% diag(svd.r$d[1:ncomp], ncomp,
                                     ncomp)
  supnorm <- apply(abs(M %*% vars), 2, max)
  q.norm <- as.numeric(quantile(supnorm, level))
  lb.norm <- mu.hat - q.norm * se
  ub.norm <- mu.hat + q.norm * se
  result <- list(x = x, y = if (keep.y) y else NULL, call = caLL,
                 model = NULL, par = NULL, nonpar = mu.hat, bandwidth = bandwidth,
                 degree = degree, level = level, teststat = NULL,
                 pnorm = NULL, pboot = NULL, qnorm = q.norm,
                 normscb = cbind(lb.norm, ub.norm), gridsize = gridsize, nrep = nrep)
  class(result) <- "SCBand"
  return(result)
}
