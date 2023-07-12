test_that("sim can be used to construct a policy data object", {

  data <- sim()

  expect_error(
    pd <- policy_data(
      type = "long",
      data = data$stage_data,
      baseline_data = data$baseline_data,
      id = "id",
      stage = "stage",
      action = "action",
      utility = "utility"
    ),
    NA
  )

  data <- sim(scale = 1/1000)

  pd <- policy_data(
    type = "long",
    data = data$stage_data,
    baseline_data = data$baseline_data,
    id = "id",
    stage = "stage",
    action = "action",
    utility = "utility"
  )

})


test_that("the utility contributions are calculated correctly", {
  K <- 15
  move_equivalent <- 700
  requisition_cost <- 14000
  in_service_premium <- 4000
  n <- 2000

  set.seed(1)
  data <- sim(K = K,
              n = n,
              move_equivalent = move_equivalent,
              requisition_cost = requisition_cost,
              in_service_premium = in_service_premium)

  pd <- policy_data(
    type = "long",
    data = data$stage_data,
    baseline_data = data$baseline_data,
    id = "id",
    stage = "stage",
    action = "action",
    utility = "utility"
  )

  #
  stage_data <- data$stage_data

  # utility contributions:
  tmp <- mean(stage_data[, .(U = sum(utility)), by = "id"]$U)
  expect_equal(tmp, mean(get_utility(pd)$U))
})

test_that("sim argument action_prob_limit", {

  data <- sim(action_prob_range = c(0.1, 0.9))

  expect_error(
    pd <- policy_data(
      type = "long",
      data = data$stage_data,
      baseline_data = data$baseline_data,
      id = "id",
      stage = "stage",
      action = "action",
      utility = "utility",
    ),
    NA
  )

})

test_that("sim handles deterministic rewards when given a Q-learning policy", {

  set.seed(1)
  data <- sim(deterministic_rewards = TRUE)

  expect_error(
    pd <- policy_data(
      type = "long",
      data = data$stage_data,
      baseline_data = data$baseline_data,
      id = "id",
      stage = "stage",
      action = "action",
      utility = "utility"
    ),
    NA
  )
  pd <- partial(pd, K = 4)

  pl <- policy_learn()
  suppressWarnings({
    po <- pl(pd, q_models = q_glm(), g_models = g_glm())
  })

  pf <- lapply(1:4, function(stage) get_policy_functions(po, stage = stage))

  set.seed(2)
  expect_error(
    suppressWarnings({
      pol_data <- sim(deterministic_rewards = TRUE, policy_functions = pf)
    }),
    NA
  )
})

