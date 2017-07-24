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
#' @param .cluster
#' @param .strat.by
#' @param .covariates
#'
#' @return
#' @export
#'
#' @examples
run_strat_reg.default <- function(.data,
                                  .formula,
                                  .cluster,
                                  .strat.by = NULL,
                                  .covariates = NULL, ...) {
  stopifnot(is.null(.strat.by) || length(.strat.by) == 1)
  stopifnot(is.null(.strat.by) || is.factor(magrittr::extract2(.data, .strat.by)))

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

  if (!is.null(.strat.by)) {
    strata.contrasts <- clean.data %>%
      select_(.dots = c(rhs.vars, .strat.by, .covariates)) %>% {
        strata.sizes <- model.matrix(as.formula(paste("~ ", .strat.by)), .) %>% colSums

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
  } else {
    design.mat <- model.matrix(.formula, clean.data)
    strata.contrasts <- NULL
    strata <- NULL
  }

  y <- model.frame(.formula, clean.data) %>% model.response() #reg.data$response
  fm <- lm.fit(design.mat, y)

  na.coef <- is.na(fm$coefficients)

  fm$coefficients %<>% magrittr::extract(!na.coef)

  fm$strat.by <- .strat.by
  fm$strata <- strata
  fm$strata.contrasts <- strata.contrasts
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
         p.value = calc.pvalue(statistic)) %>% {
    if (!is.null(fm$strat.by)) filter(., !stringr::str_detect(term, fm$strat.by)) else return(.)
  } %>%
    filter(.include_covar | !stringr::str_detect(term, "covar_")) %>%
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
    # BUGBUG Only allowing superpopulation predictions for now

    # Let's make the naming of contrasts (factors) a bit easier to parse. Changing back to default on function exit
    # old.contrasts <- getOption("contrasts")
    # old.contrasts["unordered"] <- "contr.Treatment"
    # old.options <- options(contrasts = old.contrasts)
    # on.exit(options(old.options), add = TRUE)
    #
    # newdata <- newdata %>%
    #   generate_strat_reg_data(formula(delete.response(terms(fm$formula))),
    #                           fm$strat.by,
    #                           fm$strata.constrasts,
    #                           intersect(fm$covariates, names(.))) %>%
    #   set_colnames(stringr::str_replace_all(colnames(.), setNames(c("", "(intercept)"), c(stringr::str_interp(":?${fm$strat.by}$"), "^(\\(Intercept\\))?$"))))
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

  is_mat_restrict <- !(is.character(test.list) || is.list(test.list))

  if (is_mat_restrict && ncol(test.list) != ncol(reg.output$model) - 1) {
    test.list <- matrix(0, nrow = nrow(test.list), ncol(reg.output$model) - 1) %>%
      set_colnames(names(reg.output$coefficients)) %>%
      inset(, colnames(test.list), test.list)
  }

  res <- if (!joint) {
    test.list %>% {
        if (!is_mat_restrict) {
          purrr::map(., ~ car::lht(reg.output, .x, vcov = reg.vcov, test = "F") %>% `attr<-`("linear.test", .x))
        } else {
          plyr::alply(., 1, function(test_row) car::lht(reg.output, test_row, vcov = reg.vcov, test = "F"))
        }
      } %>%
    # purrr::map(test.list, ~ car::lht(reg.output, .x, vcov = reg.vcov, test = "F") %>% `attr<-`("linear.test", .x)) %>%
    purrr::map_df(~ mutate(broom::tidy(.x),
                    linear.test = if (!is_mat_restrict) attr(.x, "linear.test") else "",
                    estimate = attr(.x, "value"),
                    std.error = sqrt(attr(.x, "vcov")))) %>%
    select(linear.test, estimate, std.error, statistic, p.value) %>%
    filter(!is.na(p.value))
  } else {
    new.class <- c("linear_test_result_joint", new.class)

    car::lht(reg.output, test.list, vcov = reg.vcov, test = "F") %>%
      broom::tidy() %>%
      filter(!is.na(p.value)) %>%
      mutate(linear.test = if (!is_mat_restrict) list(test.list) else "") %>%
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

#' Calculate the multiple hypotheses adjusted p-values
#'
#' Using the stepdown used in List, J. A., Shaikh, A. M., & Xu, Y. (2016). Multiple Hypothesis Testing in Experimental Economics (No. 21875). Code translated from the MATLAB code at https://github.com/seidelj/mht.
#'
#' @param origin_analysis_data
#' @param reg_formula
#' @param strat_by
#' @param cluster
#' @param covar
#' @param hypotheses
#' @param num_resample
#'
#' @return
#' @export
#'
#' @examples
strat_mht <- function(origin_analysis_data, reg_formula, strat_by, cluster, covar, hypotheses, num_resample = 1000, parallel = FALSE) {
  num_hypotheses <- NROW(hypotheses)

  actual_reg_res <- run_strat_reg(origin_analysis_data, reg_formula, cluster, strat_by, covar)
  actual_lht_res <- actual_reg_res %>% linear_tester(hypotheses)
  actual_ate <- actual_lht_res %>% pull(estimate)

  # if (is.character(hypotheses)) {
  # } else {
  #   hypo_col_names <- colnames(hypotheses)
  #
  #   actual_ate <- actual_reg_res %>%
  #     predict(hypotheses)
  #     # coefficients[hypo_col_names] %>%
  #     # multiply_by_matrix(hypotheses, .)
  # }

  actual_stat <- abs(actual_ate) %>% unlist

  # BUGBUG This should be done in parallel; for some reason lm.fit() gets stuck when we have something around 5000 observations.
  # lm.fit() is used by run_strat_reg().
  `%which_do%` <- if (parallel) `%dopar%` else `%do%`

  if (is.character(hypotheses)) {
    sim_stat <- foreach(resample_index = seq_len(num_resample), .combine = cbind, .errorhandling = "remove") %which_do% {
      origin_analysis_data %>% {
          if (!is.null(strat_by)) group_by_(., strat_by) else return(.)
        } %>%
        distinct_(cluster) %>%
        sample_frac(1, replace = TRUE) %>%
        ungroup() %>%
        left_join(origin_analysis_data, c(strat_by, cluster)) %>%
        run_strat_reg(reg_formula, cluster, strat_by, covar) %>%
        linear_tester(hypotheses) %>%
        select(estimate) %>%
        subtract(actual_ate) %>%
        abs() %>%
        unlist()
    } %>%
      unname()
  } else {
    hypo_col_names <- colnames(hypotheses)

    sim_stat <- foreach(resample_index = seq_len(num_resample), .combine = cbind, .errorhandling = "remove") %which_do% {
      origin_analysis_data %>% {
          if (!is.null(strat_by)) group_by_(., strat_by) else return(.)
        } %>%
        distinct_(cluster) %>%
        sample_frac(1, replace = TRUE) %>%
        ungroup() %>%
        left_join(origin_analysis_data, c(strat_by, cluster)) %>%
        run_strat_reg(reg_formula, cluster, strat_by, covar) %$%
        coefficients[hypo_col_names]
    } %>%
      magrittr::extract(, aaply(., 2, function(sim_col) all(!is.na(sim_col)))) %>% # Remove estimate cols with NAs
      multiply_by_matrix(hypotheses, .) %>%
      subtract(matrix(actual_ate, nrow = nrow(.), ncol = ncol(.))) %>%
      abs()
  }

  # Some resamples might be invalid and so we might have slightly less simulations than requested. See above foreach loop.
  num_success_sim <- ncol(sim_stat)

  actual_p <- sim_stat %>%
    is_weakly_greater_than(matrix(actual_stat, nrow(.), ncol(.))) %>%
    rowSums() %>%
    divide_by(num_success_sim) %>%
    subtract(1, .)

  sim_p <- foreach(sim_index = seq_len(num_success_sim), .combine = cbind) %do% {
    sim_stat %>%
      is_weakly_greater_than(matrix(.[, sim_index], nrow(.), ncol(.))) %>%
      rowSums() %>%
      divide_by(num_success_sim) %>%
      subtract(1, .)
  }

  p_single <- sim_p %>%
    aaply(1, . %>% sort(decreasing = TRUE)) %>%
    is_weakly_less_than(matrix(actual_p, nrow(.), ncol(.))) %>%
    alply(1, which) %>%
    map_dbl(~ if (is_empty(.x)) 1 else min(.x) / num_success_sim)

  p_single_order <- order(p_single)

  # alpha_mulm <- rep(0, num_hypotheses)

  alpha_mul <- foreach(hypo_index = 1:num_hypotheses, .combine = c) %do% {
    sim_p[p_single_order[hypo_index:num_hypotheses], , drop = FALSE] %>%
      aaply(2, max) %>%
      sort(decreasing = TRUE) %>%
      is_weakly_less_than(actual_p[p_single_order[hypo_index]]) %>%
      which() %>% {
        if (is_empty(.)) 1 else min(.) / num_success_sim
      }

    # Not yet exploiting transitivity. Code very complicated and not likely to be useful with the many TakeUp hypotheses

    # if (hypo_index == 1) {
    #   alpha_mulm[hypo_index] <- alpha_mul[hypo_index]
    # } else {
    #   for (back_hypo_index in (num_hypotheses - hypo_index + 1):1) {
    #     combn(p_single_order[hypo_index:num_hypotheses], back_hypo_index) %>%
    #       aaply(2, function(hypo_comb) {
    #         foreach(l_index = seq_len(hypo_index)) %do% {
    #           # TODO Currently assuming all hypotheses are for the same outcome and subgroup
    #
    #
    #         }
    #       })
    #   }
    # }
  }

  tibble(#hypothesis = if(is.character(hypotheses)) hypotheses else NULL,
         estimate = actual_ate,
         white_p_value = actual_lht_res %>% pull(p.value),
         unadj_p_value = p_single,
         adj_p_value = alpha_mul[order(p_single_order)]) %>%
    `attr<-`("num_resample", num_success_sim)
}
