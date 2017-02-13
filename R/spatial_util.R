#' Title
#'
#' @param in.adm
#' @param bounding.adm
#' @param max.distance
#'
#' @return
#' @export
#'
#' @examples
get.boundary.adm <- function(in.adm, bounding.adm, max.distance) {
  which.mask <- foreach::foreach(index = seq_len(nrow(in.adm)), .combine = c) %do% {
    sp::gDistance(in.adm[index, ], bounding.adm) <= max.distance
  }

  return(in.adm[which.mask, ])
}

#' Title
#'
#' @param .data
#' @param .coord.formula
#' @param .proj4string
#'
#' @return
#' @export
#'
#' @examples
convert.to.sp <- function(.data, .coord.formula, .proj4string) {
  sp::coordinates(.data) <- .coord.formula
  sp::proj4string(.data) <- .proj4string

  return(.data)
}
