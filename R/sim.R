#' @export
sim <- function(
    n = 2e3,
    K = 10,
    move_equivalent = 700,
    requisition_cost  =  14000,
    in_service_premium = 4000,
    deterministic_rewards = FALSE,
    policy_functions = NULL,
    action_prob_range = NULL,
    scale = 1
){

  data("models", package = "emrsim", envir = environment())
  list2env(models, envir = environment())

  # checks
  if(K<=1)
    stop("K must be higher than or equal to 2.")

  # extending list of policy_functions to K+1
  if (!is.null(policy_functions)){
    K_pf <- length(policy_functions)
    stopifnot(K_pf<= K)
    if(K_pf <= K){
      policy_functions <- append(policy_functions,vector(mode='list', length=(K-K_pf+1)))
    }
  }

  # baseline:
  baseline_data <- data.table(
    box = sample(model_box$action_set,
                 size = n,
                 replace = TRUE,
                 prob = model_box$tab$empir_prob)
  )
  unit <- unit_ <- rbinom(n = n,
                          size = 1,
                          prob = polle:::predict.g_empir(model_unit, baseline_data)[,1])
  unit[unit_ == 1] <- model_unit$action_set[1]
  unit[unit_ == 0] <- model_unit$action_set[2]

  baseline_data[, unit := unit]
  rm(unit_, unit)
  baseline_data[, id := 1:n]

  # stage data 1:
  stage_data_1 <- data.table(id = 1:n,
                             stage = rep(1, n))


  ## time:
  coef_model <- model_time_stage_1$coef
  rate <- 1/exp(coef_model["log(scale)"])
  shape <- exp(coef_model["log(shape)"])
  coef_model <- coef_model[!(names(coef_model) %in% c("log(scale)", "log(shape)"))]

  form <- model_time_stage_1$formula
  tt <- terms(form, data=baseline_data)
  attr(tt, "intercept") <- 1
  tt <- delete.response(tt)

  mm <- model.matrix(
    tt,
    data = baseline_data
  )
  idx_inter <- which(colnames(mm) == "(Intercept)")
  mm <- mm[,-idx_inter]

  stopifnot(
    all(colnames(mm) == names(coef_model))
  )
  mu <- c(mm %*% coef_model)
  rm(form, tt, mm, idx_inter, coef_model)

  sample_time <- (- log(runif(n)) / exp(mu))^(1/shape) / (1/rate)
  sample_age <- sample_time + 5*(365/7)
  # the first event must be before the end of the study:
  while(any(sample_age > 16.5*(365/7))){
    mu_ <- mu[sample_age > 16.5*(365/7)]
    sample_time[sample_age > 16.5*(365/7)] <-
      (- log(runif(length(mu_))) / exp(mu_))^(1/shape) / (1/rate)
    sample_age <- sample_time + 5*(365/7)
    rm(mu_)
  }
  stage_data_1[, time := sample_time]
  rm(sample_time)
  rm(rate, shape, mu)

  ## survival indicator:
  stage_data_1[, survival_indicator := FALSE]

  ## age:
  stage_data_1[, age := sample_age]
  rm(sample_age)

  ## age lag:
  stage_data_1[, age_lag := 5 * (365/7)]

  ## moves:
  lambda <- predict(model_moves,
                    newdata = cbind(baseline_data, stage_data_1),
                    type = "response")
  stage_data_1[, moves := rpois(n = n, lambda = lambda)]
  rm(lambda)

  ## cumulative moves:
  stage_data_1[, cumulative_moves := moves]

  ## cumulative moves lag:
  stage_data_1[, cumulative_moves_lag := 0]

  ## additional cost:
  stage_data_1[, cumulative_cost_lag := 0]
  sigma_ <- sigma(model_additional_cost)
  mu <- predict(model_additional_cost,
                type = "response",
                newdata = cbind(baseline_data, stage_data_1))
  sample_additional_cost <- rnorm(n = n, mean = mu, sd = sigma_)
  stage_data_1[, additional_cost := sample_additional_cost]
  rm(sample_additional_cost, sigma_, mu)

  ## work order cost:
  sigma_ <- sigma(model_work_order_cost_stage_1)
  mu <- predict(model_work_order_cost_stage_1,
                newdata = cbind(baseline_data, stage_data_1),
                type = "response")
  sample_work_order_cost <- rnorm(n = n, mean = mu, sd = sigma_)
  stage_data_1[, work_order_cost := sample_work_order_cost]
  rm(sample_work_order_cost, mu, sigma_)

  ## action lag:
  stage_data_1[, action_lag := 1]

  ## action:
  if (is.null(policy_functions)){
    prob <- predict(model_action,
                    newdata = cbind(baseline_data[, -c("id")], stage_data_1[, -c("id", "stage")]),
                    type = "response")
    if (!is.null(action_prob_range)){
      # constricts probabilities to the interval given by action_prob_range
      prob <- action_prob_range[1] * (prob < action_prob_range[1]) +
        prob * (prob >= action_prob_range[1]) * (prob <= action_prob_range[2]) +
        action_prob_range[2] * (prob > action_prob_range[2])
    }
    sample_action <- rbinom(
      n = n,
      size = 1,
      prob = prob
    )
  } else{
    H <- cbind(baseline_data[, -c("id")], stage_data_1[, -c("id", "stage")])
    if (deterministic_rewards == TRUE){
      H <- cbind(H, U_A1 = -stage_data_1[["work_order_cost"]], U_A0 = 0)
    }
    pf <- policy_functions[[1]]
    sample_action <- as.numeric(pf(H))
  }
  stage_data_1[, action := sample_action]
  rm(sample_action)

  # cumulative cost:
  stage_data_1[, cumulative_cost := action * work_order_cost + additional_cost]

  stage_data <- list()
  stage_data[[1]] <- stage_data_1
  rm(stage_data_1)
  for(k in 2:(K+1)){
    # stage data k:
    stage_data_k <- data.table(id = 1:n,
                               stage = rep(k, n))

    # age lag:
    sample_age_lag <- stage_data[[k-1]]$age
    stage_data_k[, age_lag := sample_age_lag]

    # cumulative cost lag:
    stage_data_k[
      ,cumulative_cost_lag := stage_data[[k-1]][["cumulative_cost"]]
    ]

    # cumulative moves lag:
    stage_data_k[
      ,cumulative_moves_lag := stage_data[[k-1]][["cumulative_moves"]]
    ]

    # action lag:
    sample_action_lag <- stage_data[[k-1]]$action
    stage_data_k[, action_lag := sample_action_lag]

    # work order cost lag:
    sample_work_order_cost_lag <- stage_data[[k-1]]$work_order_cost
    stage_data_k[, work_order_cost_lag := sample_work_order_cost_lag]

    ## time:
    coef_model <- model_time_stage_k$coef
    rate <- 1/exp(coef_model["log(scale)"])
    shape <- exp(coef_model["log(shape)"])
    coef_model <- coef_model[!(names(coef_model) %in% c("log(scale)", "log(shape)"))]
    form <- model_time_stage_k$formula
    environment(form) <- environment()
    tt <- terms(form, data=stage_data_k)
    attr(tt, "intercept") <- 1
    tt <- delete.response(tt)
    mm <- model.matrix(
      tt,
      data = cbind(baseline_data, stage_data_k)
    )
    idx_inter <- which(colnames(mm) == "(Intercept)")
    mm <- mm[,-idx_inter]
    stopifnot(
      all(colnames(mm) == names(coef_model))
    )
    mu <- c(mm %*% coef_model)
    rm(form, tt, mm, idx_inter, coef_model)

    sample_time <- (- log(runif(n)) / exp(mu))^(1/shape) / (1/rate)
    rm(rate, shape, mu)
    sample_time <- sample_time * sample_action_lag
    sample_age <- sample_time + sample_age_lag
    # limiting the age to 16.5*365/7 weeks
    sample_survival_indicator <- sample_age >= 16.5*365/7
    sample_age <- pmin(sample_age, 16.5*365/7)
    sample_time <- sample_age - sample_age_lag

    stage_data_k[, time := sample_time]
    stage_data_k[, survival_indicator := sample_survival_indicator]
    rm(sample_time, sample_age, sample_survival_indicator)

    ## age:
    stage_data_k[, age := time + age_lag]

    # moves:
    lambda <- predict(model_moves,
                      newdata = cbind(baseline_data, stage_data_k),
                      type = "response")
    stage_data_k[, moves := rpois(n = n, lambda = lambda) * (time>0)]
    stage_data_k[, cumulative_moves := moves + cumulative_moves_lag]
    rm(lambda)

    ## additional cost:
    sigma_ <- sigma(model_additional_cost)
    mu <- predict(model_additional_cost,
                  type = "response",
                  newdata = cbind(baseline_data, stage_data_k))

    sample_additional_cost <- rnorm(n = n, mean = mu, sd = sigma_)
    stage_data_k[, additional_cost := sample_additional_cost * (time>0)]
    rm(sigma_, mu, sample_additional_cost)

    ## work order cost:
    sigma_ <- sigma(model_work_order_cost_stage_k)
    mu <- predict(model_work_order_cost_stage_k,
                  newdata = cbind(baseline_data, stage_data_k),
                  type = "response")
    sample_work_order_cost <- rnorm(n = n, mean = mu, sd = sigma_)
    rm(sigma_,mu)
    stage_data_k[, work_order_cost := action_lag * sample_work_order_cost * (!survival_indicator)]

    ## action:
    if (is.null(policy_functions)){
      pf <- NULL
    } else{
      pf <- policy_functions[[k]]
    }
    if (is.null(pf)){
      prob <- predict(model_action,
                      newdata = cbind(baseline_data, stage_data_k),
                      type = "response")
      if (!is.null(action_prob_range)){
        # constricts probabilities to the interval given by action_prob_range
        prob <- action_prob_range[1] * (prob < action_prob_range[1]) +
          prob * (prob >= action_prob_range[1]) * (prob <= action_prob_range[2]) +
          action_prob_range[2] * (prob > action_prob_range[2])
      }
      sample_action <- rbinom(
        n = n,
        size = 1,
        prob = prob
      )
    } else {
      H <- cbind(baseline_data[, -c("id")], stage_data_k[, -c("id", "stage")])
      if (deterministic_rewards == TRUE){
        H <- cbind(H, U_A1 = -stage_data_k[["work_order_cost"]], U_A0 = 0)
      }
      sample_action <- as.numeric(pf(H))
    }
    # the K'th action is always rejected (0)
    if (k == K)
      sample_action <- sample_action * 0
    stage_data_k[, action := action_lag * sample_action * (!survival_indicator)]
    rm(sample_action)

    # cumulative cost:
    stage_data_k[, cumulative_cost := action * work_order_cost + additional_cost + cumulative_cost_lag]

    stage_data[[k]] <- stage_data_k
    rm(stage_data_k)
  }

  # combining the stage data:s
  stage_data <- rbindlist(stage_data, fill = TRUE)
  setkeyv(stage_data, c("id", "stage"))

  # removing redundant rows
  stage_data[, tmp := shift(time) == 0, by = "id"]
  stage_data <- stage_data[(is.na(tmp) | !tmp)]
  stage_data[, tmp := NULL]
  stage_data <- stage_data[!(survival_indicator == TRUE & time == 0),]

  # adding an event indicator as required by polle::policy_data
  stage_data[, event := 0]
  stage_data[time == 0, event := 1]
  stage_data[survival_indicator == TRUE, event := 1]

  # removing actions/work order cost from rows with no event
  stage_data[event == 1,  action := NA]
  stage_data[event == 1,  work_order_cost := NA]

  # sales price of reefers which are sold after the work order has been rejected.
  # reefers which are not sold receives a in service premium:
  sales_data <- stage_data[action == 0, ]
  sales_data <- merge(sales_data, baseline_data, by = "id", all.x = TRUE)
  setkeyv(sales_data, c("id", "stage"))
  mu <- predict(model_sales_price,
                newdata = sales_data,
                type = "response")
  sigma_ <- sigma(model_sales_price)
  sample_sales_price <- rnorm(n = nrow(sales_data), mean = mu, sd = sigma_)
  # bounding the sales price:
  sample_sales_price <- (sample_sales_price > 1e4) * 1e4 + (sample_sales_price <= 1e4) * sample_sales_price

  sales_data[, sales_price := sample_sales_price]
  rm(sample_sales_price)
  sales_data <- sales_data[,c("id", "sales_price")]
  sales_data[, event := 1]
  in_service_premium_data <- stage_data[survival_indicator == TRUE, ][,c("id")]
  in_service_premium_data[, sales_price := in_service_premium]
  in_service_premium_data[, event := 1]
  sales_data <- rbind(sales_data, in_service_premium_data)

  stage_data <- merge(stage_data, sales_data, all.x = TRUE, by = c("id", "event"))
  setkeyv(stage_data, c("id", "stage"))

  # calculating the utility contributions:
  stage_data[
    ,
    utility := replace(sales_price, is.na(sales_price), 0) - shift(work_order_cost, fill = 0) * action_lag + moves * move_equivalent - additional_cost,
    by = "id"
  ]
  stage_data[
    event == 1,
    utility := utility - requisition_cost
  ]

  # adding the deterministic rewards
  if (deterministic_rewards == TRUE){
    stage_data[event == 0, U_A1 := -work_order_cost]
    if (scale != 1){
      stage_data[, U_A1 := U_A1 * scale]
    }
  }

  # scaling the rewards
  if (scale != 1){
    stage_data[, utility := utility * scale]
  }

  # clean-up
  stage_data[, c(
    "age_lag",
    "action_lag",
    "work_order_cost_lag",
    "cumulative_moves_lag",
    "cumulative_cost",
    "sales_price",
    "survival_indicator"
  ) := NULL]

  return(list(
    baseline_data=baseline_data,
    stage_data=stage_data)
  )
}
