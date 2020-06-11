#' @import purrr
#' @import stats
#' @importFrom magrittr %>%
#' @importFrom utils capture.output
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))


#' @export
#' @method parallel bootstrap
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
#' @export
#' @method normal bootstrap
blbglm <- function(formula, data, m = 10, B = 5000) {
  data_list <- split_data(data, m)
  estimates <- map(
    data_list,
    ~ glm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blbglm"
  invisible(res)
}


#' split data into m parts of approximated equal sizes
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


#' compute the estimates
glm_each_subsample <- function(formula, data, n, B) {
  replicate(B, glm_each_boot(formula, data, n), simplify = FALSE)
}


#' compute the regression estimates for a blb dataset
glm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  glm1(formula, data, freqs)
}


#' estimate the regression estimates based on given the number of repetitions
glm1 <- function(formula, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  object <- glm(formula, data, weights = freqs)
  list(coef = blbcoef(object), sigma = blbsigma(object))
}


#' compute the coefficients from object
blbcoef <- function(object) {
  coef(object)
}


#' compute sigma from object
blbsigma <- function(object) {
  p <- object$rank
  y <- model.extract(object$model, "response")
  e <- fitted(object) - y
  w <- object$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' @export
#' @method print blbglm
print.blblm <- function(x, ...) {
  cat("blbglm model:", capture.output(x$formula))
  cat("\n")
}


#' @export
#' @method sigma blbglm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
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

#' @export
#' @method coef blbglm
coef.blbglm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}


#' @export
#' @method confint blbglm
confint.blbglm <- function(object, parm = NULL, level = 0.95, ...) {
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

#' @export
#' @method predict blbglm
predict.blbglm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
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
  c(object = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
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
