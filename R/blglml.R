#' @name blglml
#' @title blglml
#' @import purrr
#' @import stats
#' @importFrom magrittr %>%
#' @importFrom utils capture.output
#' @importFrom utils read.csv
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
#' @example blbglm(mpg~cyl,data = mtcars)
library(parallel)
library(future)
library(furrr)
blbglm_par <- function(formula, data, m = 10, B = 5000, Cluster) {
  data_list <- split_data(data, m)
  suppressWarnings(plan(multiprocess,workers=Cluster))
  estimates <- future_map(data_list,
                          ~glm_each_subsample(formula, data = data, n = nrow(data), B=B))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blbglm"
  invisible(res)
}
blbglm <- function(formula, data, m = 10, B = 5000) {
  data_list <- split_data(data, m)
  estimates <- map(
    data_list,
    ~ glm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blbglm"
  invisible(res)
}

#' read data in to one list
#'
#' @param filename "name"
#' @param n integer
data <- function(filename,n){
  df <- file.path(filename, list.files(filename))
  read.csv(df[n])
}



#' split data into m parts of approximated equal sizes
#'
#' @param m integer
#'
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


#' compute the estimates
#'
#' @param n integer
#' @param B integer
#'
glm_each_subsample <- function(formula, data, n, B) {
  replicate(B, glm_each_boot(formula, data, n), simplify = FALSE)
}


#' compute the regression estimates for a blb dataset
#'
#' @param n integer
#'
glm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  glm1(formula, data, freqs)
}


#' estimate the regression estimates based on given the number of repetitions
#'
#' @param freqs vector
#'
glm1 <- function(formula, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  fit <- glm(formula, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}


#' compute the coefficients from fit
#'
#' @param fit model
#'
blbcoef <- function(fit) {
  coef(fit)
}


#' compute sigma from fit
#'
#' @param fit model
#'
blbsigma <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' @export
#' @method print blbglm
print.blbglm <- function(x, ...) {
  cat("blbglm model:", capture.output(x$formula))
  cat("\n")
}

#' compute sigma from fit
#'
#' @param fit model
#'
#' @export
#' @method sigma blbglm
sigma.blbglm <- function(fit, confidence = FALSE, level = 0.95, ...) {
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
#' @param fit model
#'
#' @export
#' @method coef blbglm
coef.blbglm <- function(fit, ...) {
  est <- fit$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}

#' compute confidence interval from fit
#'
#' @param fit model
#'
#' @export
#' @method confint blbglm
confint.blbglm <- function(fit, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(fit$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- fit$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}


#' predict bootstrap from fit
#'
#' @param fit model
#'
#' @export
#' @method predict blbglm
predict.blbglm <- function(fit, new_data, confidence = FALSE, level = 0.95, ...) {
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
