#' Block bootstrap
#'
#' From https://cran.r-project.org/web/packages/broom/vignettes/bootstrapping.html
#'
#' @param df data to bootstrap
#' @param m number of samples
#' @param block string vector of IDs identifying blocks
#' @param .valid.pred predicate function to accept samples
#'
#' @return
#' @export
#'
#' @examples
block_bootstrap_ <- function(df, m, block, .valid.pred) {
  if (missing(block)) {
    n <- nrow(df)

    attr(df, "indices") <- replicate(m, sample(n, replace = TRUE) - 1, simplify = FALSE)
    attr(df, "group_sizes") <- rep(n, m)
  } else {
    block.info <- df %>%
      mutate(bootstrap.interal.index = seq_len(nrow(.)) - 1) %>%
      group_by_(.dots = block) %>%
      summarize(bn = n(),
                indices = list(bootstrap.interal.index)) %>%
      ungroup

    index.list <- replicate(m,
                            sample_n(block.info, nrow(block.info), replace = TRUE) %$% unlist(indices),
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

#' Block bootstrap
#'
#' @param df data to bootstrap
#' @param m number of samples
#' @param ... IDs identifying blocks
#' @param .valid.pred predicate function to accept samples
#'
#' @return
#' @export
#'
#' @examples
block_bootstrap <- function(df, m, ..., .valid.pred) {
  block_bootstrap_(df, m, lazyeval::lazy_dots(...), .valid.pred = .valid.pred)
}
