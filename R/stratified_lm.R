calc.pvalue <- function(t.val) {
  t.val %>%
    abs %>%
    pnorm(lower.tail = FALSE) %>%
    magrittr::multiply_by(2)
}

generate_strat_reg_data <- function(.data, .formula, .strat.by, .strata.contrasts, .covariates) {
  rhs.vars <- terms(.formula) %>%
    delete.response() %>%
    all.vars()

  stratum.contrasts.list <- setNames(list(.strata.contrasts), .strat.by)

  if (.strat.by %in% names(.data)) {
    design.mat <- update.formula(.formula, as.formula(paste("~ . *", .strat.by))) %>%
      model.matrix(.data, contrasts.arg = stratum.contrasts.list) %>%
      magrittr::extract(, stringr::str_detect(colnames(.), .strat.by)) #%>%
      # set_colnames(stringr::str_replace_all(colnames(.), sprintf("(%s)\\[T\\.([^\\]]+)\\]", paste(rhs.vars, collapse = "|")), "\\2"))
  } else {
    design.mat <- model.matrix(.formula, .data, contrasts.arg = stratum.contrasts.list)
  }

  design.mat %<>%
    set_colnames(stringr::str_replace_all(colnames(.), sprintf("(%s)\\[T\\.([^\\]]+)\\]", paste(rhs.vars, collapse = "|")), "\\2"))

  # Handling the covariates on their own because we want to demean them for each stratum. That way the intercept has a clearer
  # interpretation
  if (!purrr::is_empty(.covariates)) {
    design.mat <- .data %>%
      transmute_(.dots = c(.strat.by, setNames(.covariates, paste0("covar_", .covariates)))) %>%
      group_by_(.strat.by) %>%
      do(strat.covar.mat = model.matrix(as.formula(sprintf("~ (%s) * %s", paste(paste0("covar_", .covariates), collapse = " + "), .strat.by)),
                                        data = .,
                                        contrasts.arg = stratum.contrasts.list) %>%
          magrittr::extract(, stringr::str_detect(colnames(.), paste0("covar_.*", .strat.by))) %>%
          t %>%
          magrittr::subtract(rowMeans(.))) %>%
      ungroup %>% {
        do.call(cbind, .$strat.covar.mat)
      } %>%
      t %>%
      cbind(design.mat, .)


    # covar.row.means <- covar.design.mat %>%
    #   group_by_(.strat.by) %>%
    #   do(strat.covar.means = rowMeans(.$strat.covar.mat[[1]])) %>%
    #   ungroup
    #
    # design.mat <- covar.design.mat %>% {
    #   do.call(cbind, .$strat.covar.mat)
    #   } %>%
    #   t %>%
    #   cbind(design.mat, .)

  }

  return(design.mat)

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
  stopifnot(length(.strat.by) == 1)

  clean.data <- .data %>%
    select_(.dots = c(all.vars(.formula), .strat.by, .cluster, .covariates)) %>%
    na.omit()

  # Let's make the naming of contrasts (factors) a bit easier to parse. Changing back to default on function exit
  old.contrasts <- getOption("contrasts")
  old.contrasts["unordered"] <- "contr.Treatment"
  old.options <- options(contrasts = old.contrasts)
  on.exit(options(old.options), add = TRUE)

  rhs.vars <- terms(.formula) %>%
    delete.response() %>%
    all.vars()

  strata.contrasts <- clean.data %>%
    # select_(.dots = c(rhs.vars, "stratum", .covariates)) %>% {
    select_(.dots = c(rhs.vars, .strat.by, .covariates)) %>% {
      strata.sizes <- model.matrix(as.formula(paste("~ ", .strat.by)), .) %>% colSums

      # contr.Treatment(levels(.$stratum), contrasts = FALSE) %>%
      contr.Treatment(levels(.[[.strat.by]]), contrasts = FALSE) %>%
        magrittr::inset(1, , strata.sizes)
    }

  strata <- colnames(strata.contrasts) %>%
    stringr::str_replace("\\[(.+)\\]", "\\1")
  colnames(strata.contrasts)[1] <- ""

  strata.contrasts[1, ] %<>% magrittr::divide_by(.[1] - sum(.[-1]))
  strata.contrasts[1, -1] %<>% magrittr::multiply_by(-1)

  design.mat <- generate_strat_reg_data(clean.data, .formula, .strat.by, strata.contrasts, .covariates) %>%
    set_colnames(stringr::str_replace_all(colnames(.), setNames(c("", "(intercept)"), c(stringr::str_interp(":?${.strat.by}$"), "^$"))))

  y <- model.frame(.formula, clean.data) %>% model.response() #reg.data$response
  fm <- lm.fit(design.mat, y)

  na.coef <- is.na(fm$coefficients)

  fm$coefficients %<>% magrittr::extract(!na.coef)

  fm$strat.by <- .strat.by
  fm$strata <- strata
  fm$strata.constrasts <- strata.contrasts
  fm$cluster <- unname(unlist(clean.data[, .cluster])) #reg.data$cluster
  fm$cluster.var.name <- .cluster
  fm$model <- cbind(y, design.mat[, !na.coef])
  fm$formula <- .formula
  fm$covariates <- .covariates

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
    filter(!stringr::str_detect(term, fm$strat.by),
           .include_covar | !stringr::str_detect(term, "covar_")) %>%
    arrange(stringr::str_detect(term, "covar_")) %>% # Put the covariates last
    mutate(term = stringr::str_replace(term, "^covar_", ""))
}

#' Predict outcomes for stratified linear regressions
#'
#' @param fm
#' @param newdata
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
predict.lm_strat <- function(fm, newdata, ...) {
  if (missing(newdata)) {
    newdata <- fm$model[, -1]
  } else {
    # Let's make the naming of contrasts (factors) a bit easier to parse. Changing back to default on function exit
    old.contrasts <- getOption("contrasts")
    old.contrasts["unordered"] <- "contr.Treatment"
    old.options <- options(contrasts = old.contrasts)
    on.exit(options(old.options), add = TRUE)

    newdata <- newdata %>%
      generate_strat_reg_data(formula(delete.response(terms(fm$formula))),
                              fm$strat.by,
                              fm$strata.constrasts,
                              intersect(fm$covariates, names(.))) %>%
      set_colnames(stringr::str_replace_all(colnames(.), setNames(c("", "(intercept)"), c(stringr::str_interp(":?${fm$strat.by}$"), "^(\\(Intercept\\))?$"))))
  }


  fm$coefficients %>%
    magrittr::extract(intersect(names(.), colnames(newdata))) %>%
    magrittr::multiply_by_matrix(newdata[, names(.)], .)
}

#' Calculate the R Squared for stratified regressions
#'
#' @param fm stratified regression results
#' @param adjusted
#'
#' @return
#' @export
#'
#' @examples
r_squared <- function(fm, adjusted = TRUE) {
  # Taken from summary.lm()
  mss <- sum((fm$fitted.values - mean(fm$fitted.values))^2)
  rss <- sum(fm$residuals^2)

  r.squared <- mss / (mss + rss)

  if (adjusted) {
    r.squared <- 1 - (1 - r.squared) * ((length(fm$fitted.values) - 1)/fm$df.residual)
  }

  return(r.squared)
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

  reg.vcov <- vcov_clx(reg.output)

  res <- if (!joint) {
    purrr::map(test.list, ~ car::lht(reg.output, .x, vcov = reg.vcov, test = "F") %>% `attr<-`("linear.test", .x)) %>%
    purrr::map_df(~ mutate(broom::tidy(.x),
                    linear.test = attr(.x, "linear.test"),
                    estimate = attr(.x, "value"),
                    std.error = sqrt(attr(.x, "vcov")))) %>%
    select(linear.test, estimate, std.error, statistic, p.value) %>%
    filter(!is.na(p.value))
  } else {
    new.class <- c("linear_test_result_joint", new.class)

    car::lht(reg.output, test.list, vcov = reg.vcov, test = "F") %>%
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
