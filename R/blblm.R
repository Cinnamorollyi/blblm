#' @import purrr
#' @import stats
#' @importFrom utils capture.output
#' @importFrom utils read.csv
#' @importFrom magrittr %>%
#'
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))


#' Using parallel/normal function to do bootstrap
#'
#' @param data file
#' @param Cluster cores
#'
#' @return list
#' @export
#' @method two bootstrap
#' @example blblm(mpg~cyl,data = mtcars)
library(parallel)
library(future)
library(furrr)
blblm_par <- function(formula, data, m = 10, B = 5000, Cluster) {
  data_list <- split_data(data, m)
  suppressWarnings(plan(multiprocess,workers=Cluster))
  estimates <- future_map(data_list,
                          ~lm_each_subsample(formula, data = data, n = nrow(data), B=B))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}

blblm <- function(formula, data, m = 10, B = 5000) {
  data_list <- split_data(data, m)
  estimates <- map(
    data_list,
    ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}

#' read data in to one list
data <- function(filename,n){
  df <- file.path(filename, list.files(filename))
  read.csv(df[n])
}


#' split data into m parts of approximated equal sizes
split_data <- function(data, m) {
  df <- vroom(list_of_files, .id = "FileName")
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


#' compute the estimates
lm_each_subsample <- function(formula, data, n, B) {
  replicate(B, lm_each_boot(formula, data, n), simplify = FALSE)
}


#' compute the regression estimates for a blb dataset
lm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  lm1(formula, data, freqs)
}


#' estimate the regression estimates based on given the number of repetitions
lm1 <- function(formula, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  fit <- lm(formula, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}


#' compute the coefficients from fit
blbcoef <- function(fit) {
  coef(fit)
}


#' compute sigma from fit
blbsigma <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}

#' @export
#' @return character
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}


#' compute sigma from fit
#'
#' @param object model
#' @param confidence boolean vector
#' @param level confidencial level
#' @param ... est
#' @return double
#' @export
#' @method sigma blbglm
sigma.blblm <- function(fit, confidence = FALSE, level = 0.95, ...) {
  est <- fit$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' compute coef from fit
#'
#' @param object model
#' @param ... est
#' @return tibble
#' @export
#' @method coef blbglm
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}


#' compute confidence interval from fit
#' @param object model
#' @param parm boolean vector
#' @param level cofidential level
#' @param ... est
#' @return out tibble
#' @export
#' @method confint blbglm
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' compute confidence interval from fit
#'
#' @param fit model
#' @param new_data matrics/vector
#' @param confidence boolean vector
#' @param level confidential level
#' @param ... est
#' @return tibble
#' @export
#' @method confint blbglm
predict.blblm <- function(fit, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- fit$estimates
  X <- model.matrix(reformulate(attr(terms(fit$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
               apply(1, mean_lwr_upr, level = level) %>%
               t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}