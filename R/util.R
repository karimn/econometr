#' Title
#'
#' @param .values
#'
#' @return
#' @export
#'
#' @examples
is_outlier <- function(.values) {
  .values %>% { . < quantile(., 0.25) - 1.5 * IQR(.) | . > quantile(., 0.75) + 1.5 * IQR(.) }
}
