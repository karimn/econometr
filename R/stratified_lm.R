calc.pvalue <- function(t.val) {
  t.val %>%
    abs %>%
    pnorm(lower.tail = FALSE) %>%
    magrittr::multiply_by(2)
}

#' @export
run_strat_reg  <- function(.data, ...) UseMethod("run_strat_reg")

#' Run Stratified Regression
#'
#' @param .data
#' @param .formula
#' @param .strat.by
#' @param .cluster
#' @param .covariates
#'
#' @return
#' @export
#'
#' @examples
run_strat_reg.default <- function(.data,
                                  .formula,
                                  .strat.by,
                                  .cluster,
                                  .covariates = NULL, ...) {
  clean.data <- .data %>%
    select_(.dots = c(union(c(all.vars(.formula), .strat.by), .covariates), .cluster)) %>%
    na.omit %>%
    tidyr::unite_("stratum", from = .strat.by, sep = ".", remove = TRUE) %>%
    mutate(stratum = factor(stratum))

  get.strat.design.matrix <- function() {
    strata.contrasts <- clean.data %>%
    select_(.dots = c(all.vars(.formula)[-1], "stratum", .covariates)) %>% {
      strata.sizes <- model.matrix(~ stratum, .) %>% colSums
      contrasts(.$stratum, contrasts = FALSE) %>%
        magrittr::inset(1, , strata.sizes)
    }

    colnames(strata.contrasts)[1] <- ""

    strata.contrasts[1, ] %<>% magrittr::divide_by(.[1] - sum(.[-1]))
    strata.contrasts[1, -1] %<>% magrittr::multiply_by(-1)

    strat.design.mat <- clean.data %>%
      model.matrix(~ stratum, ., contrasts.arg = list(stratum = strata.contrasts)) %>%
      magrittr::extract(, -1)

    treat.design.mat <- model.matrix(.formula, clean.data) %>%
      magrittr::extract(, -1) %>%
      set_colnames(stringr::str_replace_all(colnames(.), paste(all.vars(.formula)[-1], collapse = "|"), "")) %>% {
        new.col.names <- paste(rep(colnames(.), each = ncol(strat.design.mat)), colnames(strat.design.mat), sep = "_")
        plyr::alply(., 2, function(.treat.col) strat.design.mat * .treat.col) %>% { do.call(cbind, .) } %>%
          set_colnames(new.col.names)
      }


    if (!is.null(.covariates)) {
      covar.design.mat <- clean.data %>%
        group_by(stratum) %>%
        do(strat.covar.mat = {
          model.matrix(as.formula(sprintf("~ %s", paste(.covariates, collapse = " + "))), .) %>%
            magrittr::extract(, -1) %>%
            t %>%
            magrittr::subtract(rowMeans(.))
        }) %>%
        ungroup %>% {
          do.call(cbind, .$strat.covar.mat)
        } %>%
        t %>%
        set_colnames(paste0("covar_", colnames(.))) %>% {
          new.col.names <- paste(rep(colnames(.), each = ncol(strat.design.mat)), colnames(strat.design.mat), sep = "_")

          unname(.) %>%
            plyr::alply(2, function(.covar.col) strat.design.mat * .covar.col) %>% { do.call(cbind, .) } %>%
            set_colnames(new.col.names)
        }
      design.mat <- cbind(strat.design.mat, treat.design.mat, covar.design.mat)
    } else {
      design.mat <- cbind(strat.design.mat, treat.design.mat)
    }

    design.mat %>%
      set_colnames(stringr::str_replace_all(colnames(.), list("_?stratum$" = "", "^$" = "(intercept)")))
  }

  design.mat <- get.strat.design.matrix()

  y <- model.frame(.formula, clean.data) %>% model.response()

  fm <- lm.fit(design.mat, y)

  na.coef <- is.na(fm$coefficients)

  fm$coefficients %<>% magrittr::extract(!na.coef)

  fm$cluster <- unname(unlist(clean.data[, .cluster]))
  fm$model <- cbind(y, design.mat[, !na.coef])

  class(fm) <- "lm_strat"

  return(fm)
}

#' Title
#'
#' @param fm
#'
#' @return
#' @export
#'
#' @examples
estfun.lm_strat <- function(fm) {
  fm$model %>%
    magrittr::extract(, colnames(.) != "y") %>%
    magrittr::multiply_by(fm$residuals)
}

#' Title
#'
#' @param fm
#'
#' @return
#' @export
#'
#' @examples
bread.lm_strat <- function(fm) {
  fm$model %>%
    magrittr::extract(, colnames(.) != "y") %>%
    crossprod %>%
    solve %>%
    magrittr::multiply_by(nrow(fm$model))
}

#' Title
#'
#' @param fm
#'
#' @return
#' @export
#'
#' @examples
vcov_clx.lm_strat <- function(fm, ...) {
  vcov_clx.default(fm, fm$cluster)
}

#' Tidy stratified regression results
#'
#' @param fm
#' @param .include_covar
#'
#' @return
#' @export
#'
#' @examples
tidy.lm_strat <- function(fm, .include_covar = FALSE, ...) {
  tibble(term = names(fm$coefficients),
         estimate = fm$coefficients,
         std.error = vcov_clx(fm) %>% diag %>% sqrt,
         statistic = fm$coefficients / std.error,
         p.value = calc.pvalue(statistic)) %>%
    filter(!stringr::str_detect(term, "stratum"),
           .include_covar | !stringr::str_detect(term, "covar_")) %>%
    arrange(stringr::str_detect(term, "covar_")) %>% # Put the covariates last
    mutate(term = stringr::str_replace(term, "^covar_", ""))
}

#' Title
#'
#' @param reg.output
#' @param test.list
#' @param joint
#'
#' @return
#' @export
#'
#' @examples
linear_tester <- function(reg.output, test.list, joint = FALSE) {
  new.class <- "linear_test_result"

  res <- if (!joint) {
    purrr::map(test.list, ~ car::lht(reg.output, .x, vcov = vcov_clx, test = "F") %>% `attr<-`("linear.test", .x)) %>%
    purrr::map_df(~ mutate(broom::tidy(.x),
                    linear.test = attr(.x, "linear.test"),
                    estimate = attr(.x, "value"),
                    std.error = sqrt(attr(.x, "vcov")))) %>%
    select(linear.test, estimate, std.error, statistic, p.value) %>%
    filter(!is.na(p.value))
  } else {
    new.class <- c("linear_test_result_joint", new.class)

    car::lht(reg.output, test.list, vcov = vcov_clx, test = "F") %>%
      broom::tidy() %>%
      filter(!is.na(p.value)) %>%
      mutate(linear.test = list(test.list)) %>%
      select(linear.test, statistic, p.value)
  }

  class(res) <- c(new.class, oldClass(res))

  return(res)
}

#' Title
#'
#' @param .res
#'
#' @return
#' @export
#'
#' @examples
tidy.linear_test_result_joint <- function(.res) {
  .res %>%
    rowwise %>%
    mutate(linear.test = sprintf("(%s)", paste(linear.test, collapse = ") & ("))) %>%
    ungroup
}
