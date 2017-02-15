#' Identify outlier values
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

#' Make data wide for multiple columns
#'
#' @param .data
#' @param .grp.id
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
multi_spread <- function(.data, .grp.id, ...) {
  .data %>%
    tidyr::gather_("multi_spread_temp_key", "multi_spread_temp_val",
                   gather_cols = unname(dplyr::select_vars(colnames(.), ...))) %>%
    tidyr::unite_("multi_spread_temp", c("multi_spread_temp_key", .grp.id)) %>%
    tidyr::spread(multi_spread_temp, multi_spread_temp_val)
}
