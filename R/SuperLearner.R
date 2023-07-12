#' @export
unlock_environment <- function(env){
  return (new.env(parent=env))
}

#' @export
strip_gam <- function(op){
  op$y <- NULL
  op$model <- NULL
  op$fitted.values <- NULL
  attr(op$terms, ".Environment") <- NULL
  attr(op$formula, ".Environment") <- NULL

  if (!is.null(op$dinfo)){
    environment(op$dinfo$gp$pred.formula) <- NULL
    environment(op$dinfo$gp$fake.formula) <- NULL
    environment(op$dinfo$gp$pf) <- NULL
  }

  return(op)
}

#' @export
strip_glm <- function(cm){
  cm$y <- NULL
  cm$model <- NULL

  cm$residuals <- NULL
  cm$fitted.values <- NULL
  cm$effects <- NULL
  cm$qr$qr <- NULL
  cm$linear.predictors <- NULL
  cm$weights <- NULL
  cm$prior.weights <- NULL
  cm$data <- NULL

  attr(cm$terms, ".Environment") <- NULL
  attr(cm$formula, ".Environment") <- NULL

  return(cm)
}

#' @export
SL.glm_q <- function (Y, X, newX, family, obsWeights, model = FALSE, ...)
{

  if (is.matrix(X)) {
    X = as.data.frame(X)
  }

  form <- Y ~ A1 * .

  fit.glm <- suppressWarnings(glm(form, data = X, family = family, weights = obsWeights,
                                  model = model))

  if (is.matrix(newX)) {
    newX <- as.data.frame(newX)
  }
  pred <- predict(fit.glm, newdata = newX, type = "response")

  # removing the environment from the object
  fit.glm <- strip_glm(fit.glm)

  fit <- list(object = fit.glm)
  class(fit) <- "SL.glm"
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' @export
SL.glm_q_unit <- function (Y, X, newX, family, obsWeights, model = FALSE, ...)
{

  if (is.matrix(X)) {
    X = as.data.frame(X)
  }

  form <- Y ~ A1 * unit2 * .

  fit.glm <- suppressWarnings(glm(form, data = X, family = family, weights = obsWeights,
                                  model = model))

  if (is.matrix(newX)) {
    newX <- as.data.frame(newX)
  }
  pred <- predict(fit.glm, newdata = newX, type = "response")

  # removing the environment from the object
  fit.glm <- strip_glm(fit.glm)

  fit <- list(object = fit.glm)
  class(fit) <- "SL.glm"
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' @export
SL.glm_form <- function (Y, X, newX, family, obsWeights, model = FALSE, formula = ~., ...){
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
  tt <- terms(formula, data = X)
  form <- reformulate(attr(tt, "term.labels"), response = "Y")

  fit.glm <- suppressWarnings(glm(formula = form,
                                  data = X,
                                  family = family,
                                  weights = obsWeights,
                                  model = model))

  if (is.matrix(newX)) {
    newX <- as.data.frame(newX)
  }
  pred <- predict(fit.glm, newdata = newX, type = "response")

  # removing the environment from the object
  fit.glm <- strip_glm(fit.glm)

  fit <- list(object = fit.glm)
  class(fit) <- "SL.glm_form"
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' @export
predict.SL.glm_form <- function (object, newdata, ...)
{
  if (is.matrix(newdata)) {
    newdata = as.data.frame(newdata)
  }

  fit <- object$object

  # adding the splines package to the model terms environment
  attr(fit$terms, ".Environment") <- unlock_environment(as.environment("package:splines"))

  pred <- predict(object = fit, newdata = newdata,
                  type = "response")
  pred
}

#' @export
SL.glmnet_form <- function (Y, X, newX, family, obsWeights, id, alpha = 1, nfolds = 10,
                            formula, nlambda = 100, useMin = TRUE, loss = "deviance", ...){
  if (!require("glmnet")) {
    stop("SL.glmnet_form requires the glmnet package, but it isn't available")
  }
  if (!require("splines")) {
    stop("SL.glmnet_form requires the splines package, but it isn't available")
  }

  tt <- terms(formula, data = X, simplify = TRUE)
  if (all(newX$cumulative_cost_lag == 0)){
    tt <- drop.terms(tt,
                     which(grepl("cumulative_cost_lag", attr(tt,"term.labels"))))
  }
  form <- reformulate(attr(tt, "term.labels"))

  des <- polle:::get_design(form, data=X)
  fitCV <- glmnet::cv.glmnet(x = des$x, y = Y, weights = obsWeights,
                             lambda = NULL, type.measure = loss, nfolds = nfolds,
                             family = family$family, alpha = alpha, nlambda = nlambda,
                             ...)

  newx <- polle:::apply_design(des, newX)
  pred <- predict(fitCV, newx = newx, type = "response", s = ifelse(useMin,
                                                                    "lambda.min", "lambda.1se"))
  des$x <- NULL
  attr(des$terms, ".Environment") <- NULL

  fit <- list(object = fitCV, useMin = useMin, design = des)
  class(fit) <- "SL.glmnet_q"
  out <- list(pred = pred, fit = fit)
  return(out)
}



#' @export
SL.glmnet_q <- function (Y, X, newX, family, obsWeights, id, alpha = 1, nfolds = 10,
                         nlambda = 100, useMin = TRUE, loss = "deviance", ...){

  if (!require("glmnet")) {
    stop("SL.glmnet_q requires the glmnet package, but it isn't available")
  }
  if (!require("splines")) {
    stop("SL.glmnet_q requires the splines package, but it isn't available")
  }

  if (!all(newX$cumulative_cost_lag == 0)){ # stage 1 has no cumulative cost lag
    form <- ~ A1 * (
      ns(age, df = 5)+
        ns(work_order_cost, df = 5) +
        ns(cumulative_moves, df = 5) +
        ns(cumulative_cost_lag, df = 5) +
        age:work_order_cost +
        age:cumulative_cost_lag +
        .
    )
  } else {
    form <- ~ A1 * (
      ns(age, df = 5)+
        ns(work_order_cost,df = 5) +
        ns(cumulative_moves, df = 5) +
        age:work_order_cost +
        .
    )
  }

  des <- polle:::get_design(form, data=X)
  fitCV <- glmnet::cv.glmnet(x = des$x, y = Y, weights = obsWeights,
                             lambda = NULL, type.measure = loss, nfolds = nfolds,
                             family = family$family, alpha = alpha, nlambda = nlambda,
                             ...)

  newx <- polle:::apply_design(des, newX)
  pred <- predict(fitCV, newx = newx, type = "response", s = ifelse(useMin,
                                                                    "lambda.min", "lambda.1se"))
  des$x <- NULL
  attr(des$terms, ".Environment") <- NULL

  fit <- list(object = fitCV, useMin = useMin, design = des)
  class(fit) <- "SL.glmnet_q"
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' @export
SL.glmnet_q_unit <- function (Y, X, newX, family, obsWeights, id, alpha = 1, nfolds = 10,
                         nlambda = 100, useMin = TRUE, loss = "deviance", ...){

  if (!require("glmnet")) {
    stop("SL.glmnet_q_ns requires the glmnet package, but it isn't available")
  }
  if (!require("splines")) {
    stop("SL.glmnet_q_ns requires the splines package, but it isn't available")
  }

  if (!all(newX$cumulative_cost_lag == 0)){ # stage 1 has no cumulative cost lag
    form <- ~ A1 * unit2 * (
      ns(age, df = 5)+
        ns(work_order_cost, df = 5) +
        ns(cumulative_moves, df = 5) +
        ns(cumulative_cost_lag, df = 5) +
        age:work_order_cost +
        age:cumulative_cost_lag +
        .
    )
  } else {
    form <- ~ A1 * unit2 * (
      ns(age, df = 5)+
        ns(work_order_cost, df = 5) +
        ns(cumulative_moves, df = 5) +
        age:work_order_cost +
        .
    )
  }

  des <- polle:::get_design(form, data=X)
  fitCV <- glmnet::cv.glmnet(x = des$x, y = Y, weights = obsWeights,
                             lambda = NULL, type.measure = loss, nfolds = nfolds,
                             family = family$family, alpha = alpha, nlambda = nlambda,
                             ...)

  newx <- polle:::apply_design(des, newX)
  pred <- predict(fitCV, newx = newx, type = "response", s = ifelse(useMin,
                                                                    "lambda.min", "lambda.1se"))
  des$x <- NULL
  attr(des$terms, ".Environment") <- NULL

  fit <- list(object = fitCV, useMin = useMin, design = des)
  class(fit) <- "SL.glmnet_q"
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' @export
predict.SL.glmnet_q <- function (object, newdata, remove_extra_cols = T, add_missing_cols = T,
                                 ...)
{
  if (!require("glmnet")) {
    stop("SL.glmnet_q requires the glmnet package, but it isn't available")
  }
  if (!require("splines")) {
    stop("SL.glmnet_q requires the splines package, but it isn't available")
  }

  design <- object$design
  # adding the splines package to the terms environment
  attr(design$terms, ".Environment") <- unlock_environment(as.environment("package:splines"))
  newx <- polle:::apply_design(design, newdata)
  pred <- predict(object$object, newx = newx, type = "response",
                  s = ifelse(object$useMin, "lambda.min", "lambda.1se"))
  return(pred)
}

#' @export
SL.gam_mgcv_form <- function (Y, X, newX, family, obsWeights, formula=Y~., discrete = FALSE,...){

  if (!require("mgcv")) {
    stop("SL.gam_mgcv_form requires the mgcv package, but it isn't available")
  }
  tt <- terms(formula, data = X, simplify = TRUE)
  form <- reformulate(attr(tt, "term.labels"), response = "Y")

  fit.gam <- mgcv::bam(
    formula = form,
    data = X,
    family = family,
    weights = obsWeights,
    discrete = discrete
  )
  fit.gam <- strip_gam(fit.gam)

  pred <- mgcv::predict.gam(fit.gam, newdata = newX, type = "response")

  fit <- list(object = fit.gam)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.gam_mgcv")
  return(out)
}

#' @export
SL.gam_mgcv_q <- function (Y, X, newX, family, obsWeights, ...){

  if (!require("mgcv")) {
    stop("SL.glmnet_q_ns requires the glmnet package, but it isn't available")
  }

  if (!all(newX$cumulative_cost_lag == 0)){ # stage 1 has no cumulative cost lag
    form <- Y ~
      s(age, by = factor(A1, level = c("0", "1"))) +
      s(work_order_cost, by = factor(A1, level = c("0", "1"))) +
      s(cumulative_moves, by = factor(A1, level = c("0", "1"))) +
      s(cumulative_cost_lag, by = factor(A1, level = c("0", "1"))) +
      A1 * (
        age:work_order_cost +
          age:cumulative_cost_lag +
          region2 + region3 + region4 + region5 +
          additional_cost +
          box2 + unit2
      )
  } else {
    form <- Y ~
      s(age, by = factor(A1, level = c("0", "1"))) +
      s(work_order_cost, by = factor(A1, level = c("0", "1"))) +
      s(cumulative_moves, by = factor(A1, level = c("0", "1"))) +
      A1 * (
        age:work_order_cost +
          region2 + region3 + region4 + region5 +
          additional_cost +
          box2 + unit2
      )
  }

  fit.gam <- mgcv::bam(form, data = X, family = family,
                       weights = obsWeights)

  # stripping model
  fit.gam <- strip_gam(fit.gam)

  pred <- mgcv::predict.gam(fit.gam, newdata = newX, type = "response")

  fit <- list(object = fit.gam)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.gam_mgcv")
  return(out)
}

#' @export
SL.gam_mgcv_q_2 <- function (Y, X, newX, family, obsWeights, ...){

  if (!require("mgcv")) {
    stop("SL.glmnet_q_ns requires the glmnet package, but it isn't available")
  }

  if (!all(newX$cumulative_cost_lag == 0)){ # stage 1 has no cumulative cost lag
    form <- Y ~
      s(age, by = factor(A1, level = c("0", "1"))) +
      s(work_order_cost, by = factor(A1, level = c("0", "1"))) +
      s(cumulative_moves, by = factor(A1, level = c("0", "1"))) +
      s(cumulative_cost_lag, by = factor(A1, level = c("0", "1"))) +
      A1 * (
        age:work_order_cost +
          age:cumulative_cost_lag +
          additional_cost +
          box2 + unit2
      )
  } else {
    form <- Y ~
      s(age, by = factor(A1, level = c("0", "1"))) +
      s(work_order_cost, by = factor(A1, level = c("0", "1"))) +
      s(cumulative_moves, by = factor(A1, level = c("0", "1"))) +
      A1 * (
        age:work_order_cost +
          additional_cost +
          box2 + unit2
      )
  }

  fit.gam <- mgcv::bam(form, data = X, family = family,
                       weights = obsWeights)
  # stripping model
  fit.gam <- strip_gam(fit.gam)

  pred <- mgcv::predict.gam(fit.gam, newdata = newX, type = "response")

  fit <- list(object = fit.gam)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.gam_mgcv")
  return(out)
}

#' @export
SL.gam_mgcv_q_unit <- function (Y, X, newX, family, obsWeights, ...){

  if (!require("mgcv")) {
    stop("SL.glmnet_q_ns requires the glmnet package, but it isn't available")
  }

  if (!all(newX$cumulative_cost_lag == 0)){ # stage 1 has no cumulative cost lag
    form <- Y ~
      s(age, by = factor(paste(A1, unit2, sep = ""))) +
      s(work_order_cost, by = factor(paste(A1, unit2, sep = ""))) +
      s(cumulative_moves, by = factor(paste(A1, unit2, sep = ""))) +
      s(cumulative_cost_lag, by = factor(paste(A1, unit2, sep = ""))) +
      A1 * box2
  } else {
    form <- Y ~
      s(age, by = factor(paste(A1, unit2, sep = ""))) +
      s(work_order_cost, by = factor(paste(A1, unit2, sep = ""))) +
      s(cumulative_moves, by = factor(paste(A1, unit2, sep = ""))) +
      A1 * box2
  }

  fit.gam <- mgcv::bam(
    form,
    data = X,
    family = family,
    weights = obsWeights,
    discrete = TRUE
  )
  # stripping model
  fit.gam <- strip_gam(fit.gam)

  pred <- mgcv::predict.gam(fit.gam, newdata = newX, type = "response")

  fit <- list(object = fit.gam)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.gam_mgcv")
  return(out)
}

#' @export
SL.gam_mgcv_q_unit_te <- function (Y, X, newX, family, obsWeights, ...){

  if (!require("mgcv")) {
    stop("SL.glmnet_q_ns requires the glmnet package, but it isn't available")
  }

  if (!all(newX$cumulative_cost_lag == 0)){ # stage 1 has no cumulative cost lag
    form <- Y ~
      te(age, work_order_cost, by = factor(paste(A1, unit2, sep = ""))) +
      ti(age, I(cumulative_cost_lag + additional_cost), by = factor(paste(A1, unit2, sep = ""))) +
      s(I(cumulative_cost_lag + additional_cost), by = factor(paste(A1, unit2, sep = ""))) +
      s(cumulative_moves, by = factor(paste(A1, unit2, sep = ""))) +
      A1 * box2
  } else {
    form <- Y ~
      te(age, work_order_cost, by = factor(paste(A1, unit2, sep = ""))) +
      ti(age, additional_cost, by = factor(paste(A1, unit2, sep = ""))) +
      s(additional_cost, by = factor(paste(A1, unit2, sep = ""))) +
      s(cumulative_moves, by = factor(paste(A1, unit2, sep = ""))) +
      A1 * box2
  }

  fit.gam <- mgcv::bam(
    form,
    data = X,
    family = family,
    weights = obsWeights,
    discrete = TRUE
  )
  # stripping model
  fit.gam <- strip_gam(fit.gam)

  pred <- mgcv::predict.gam(fit.gam, newdata = newX, type = "response")

  fit <- list(object = fit.gam)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.gam_mgcv")
  return(out)
}


#' @export
predict.SL.gam_mgcv <- function (object, newdata, ...)
{
  if(!require("mgcv"))
    Stop("SL.gam_mgcv requires the mgcv package, but it isn't available")

  pred <- mgcv::predict.gam(object = object$object, newdata = newdata,
                            type = "response")

  return(pred)
}

#' @export
SL.gam_mgcv_g <- function (Y, X, newX, family, obsWeights, k = 3, bs = "cr", maxit = 100, ...)
{

  if (!require("mgcv")) {
    stop("SL.gam.mgcv requires the mgcv package, but it isn't available")
  }
  gam.model <- Y ~ te(age, work_order_cost, k = 3) + box2 + unit2 + region2 + region3 + region4 + region5

  fit.gam <- mgcv::bam(gam.model, data = X, family = family,
                       control = mgcv::gam.control(maxit = maxit),
                       weights = obsWeights)
  # stripping model
  fit.gam <- strip_gam(fit.gam)

  pred <- mgcv::predict.gam(fit.gam, newdata = newX, type = "response")

  fit <- list(object = fit.gam)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.gam_mgcv")
  return(out)
}

#' @export
SL.gam_mgcv_g_2 <- function (Y, X, newX, family, obsWeights, k = NA, bs = "cr", maxit = 200, ...)
{

  if (!require("mgcv")) {
    stop("SL.gam.mgcv requires the mgcv package, but it isn't available")
  }
  gam.model <- Y ~ te(age, work_order_cost)

  fit.gam <- mgcv::bam(gam.model, data = X, family = family,
                       control = mgcv::gam.control(maxit = maxit),
                       weights = obsWeights)
  # stripping model
  fit.gam <- strip_gam(fit.gam)

  pred <- mgcv::predict.gam(fit.gam, newdata = newX, type = "response")

  fit <- list(object = fit.gam)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.gam_mgcv")
  return(out)
}

#' @export
predict.SL.gam_mgcv <- function (object, newdata, ...)
{
  if(!require("mgcv"))
    Stop("SL.gam_mgcv requires the mgcv package, but it isn't available")

  pred <- mgcv::predict.gam(object = object$object, newdata = newdata,
                            type = "response")

  return(pred)
}
