test_that("onerun_policy_eval",{

  expect_no_error(
    suppressWarnings({
      ope <- onerun_policy_eval(
        sim = sim,
        M = 1,
        type = "or",
        n = 2e3,
        n_approx = 1e4,
        K = 3,
        K_partial = 2,
        policy_learn = policy_learn(),
        q_models = q_glm(),
        save_q_functions = FALSE
      )
    })

  )
  expect_true(is.na(ope["value_estimate_vcov"]))

  expect_no_error(
    suppressWarnings({
      ope <- onerun_policy_eval(
        sim = sim,
        M = 1,
        type = "dr",
        n = 2e3,
        n_approx = 1e4,
        K = 3,
        K_partial = 2,
        policy_learn = policy_learn(),
        q_models = q_glm(),
        g_models = g_glm(),
        save_q_functions = FALSE
      )
    })
  )
  expect_true(!is.na(ope["value_estimate_vcov"]))

})
