#' Block bootstrap
#'
#' From https://cran.r-project.org/web/packages/broom/vignettes/bootstrapping.html
#'
#' @param df
#' @param m
#' @param block
#' @param .valid.pred
#'
#' @return
#' @export
#'
#' @examples
block_bootstrap <- function(df, m, block, .valid.pred) {
  if (missing(block)) {
    n <- nrow(df)

    attr(df, "indices") <- replicate(m, sample(n, replace = TRUE) - 1, simplify = FALSE)
    attr(df, "group_sizes") <- rep(n, m)
  } else {
    block.info <- df %>%
      mutate(bootstrap.interal.index = seq_len(nrow(.)) - 1) %>%
      group_by_(block) %>%
      summarize(bn = n(),
                indices = list(bootstrap.interal.index))

    index.list <- replicate(m,
                            sample_n(block.info, nrow(block.info), replace = TRUE) %$% indices %>% unlist,
                            simplify = FALSE)

    if (!missing(.valid.pred)) {
      index.list %<>% purrr::keep(~.valid.pred(df[. + 1, ]))

      m <- length(index.list)
    }

    attr(df, "indices") <- index.list
    attr(df, "group_sizes") <- purrr::map_int(index.list, length)
  }

  attr(df, "drop") <- TRUE
  attr(df, "biggest_group_size") <- n
  attr(df, "labels") <- data.frame(replicate = 1:m)
  attr(df, "vars") <- list(quote(replicate))
  class(df) <- c("grouped_df", "tbl_df", "tbl", "data.frame")

  return(df)
}
