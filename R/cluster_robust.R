#' @export
vcov_clx <- function(fm, cluster, ...) UseMethod("vcov_clx")

#' Calculate the Cluster Robust Variance Covariance Matrix
#'
#' Taken from ivpack::clx()
#'
#' @param fm
#' @param cluster
#'
#' @return
#' @export
#'
#' @examples
vcov_clx.default <- function(fm, cluster, ...) {
  dfcw <- 1
  M <- NROW(unique(cluster))
  N <- NROW(cluster)

  if (N != nrow(fm$model)) { # Need to match cluster to fitted model data
    cluster %<>%
      merge(fm$model[, NULL], by="row.names", all=FALSE) %>%
      select(-Row.names)

    N <- NROW(cluster)

    if (N != nrow(fm$model)) {
      stop("Cannot match clusters with fitted model's data")
    }
  }

  dfc <- (M/(M - 1)) * ((N - 1)/(N - fm$rank))
  u <- apply(sandwich::estfun(fm), 2, . %>% tapply(cluster, sum))

  dfc * sandwich::sandwich(fm, meat. = crossprod(u)/N) * dfcw
}
