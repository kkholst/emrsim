name <- "pe_blip_(q_sl)_(g_sl)_(blip_sl_glm)"

library("polle")
library("emrsim")
library("future.apply")
library("progressr")

set.seed(423)
n0 <- 2e4
n_approx0 <- 1e5
K0 <- 16
K_partial0 <- 8
M0 <- 1
R0 <- 100
workers <- 20
alpha0 <- c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1)

# Q-models ----------------------------------------------------------------

## xgb
xgboost_tune_q <- list(max_depth = c(3,4),
                       shrinkage = c(0.05,0.1),
                       ntrees  = c(250, 300))
xgboost_learners_q = create.Learner(
  "SL.xgboost",
  tune = xgboost_tune_q,
  detailed_names = TRUE,
  name_prefix = "xgb"
)

q_models <- q_sl(
  SL.library = c(
    "SL.mean",
    "SL.glm_q",
    "SL.glmnet_q_unit",
    "SL.gam_mgcv_q_unit",
    "SL.gam_mgcv_q_unit_te",
    "xgb_4_0.05_300"
  ),
  cvControl = SuperLearner.CV.control(V = 4L),
  env = environment()
)

# g-models ----------------------------------------------------------------

g_models <- g_sl(
  SL.library = c(
    "SL.mean",
    "SL.glm",
    "SL.gam_mgcv_g_2"
  ),
  cvControl = SuperLearner.CV.control(V = 2L),
  env = parent.frame()
)

# blip-models ----------------------------------------------------------------

blip_formulas <- list(
  form1 = ~ .,
  form2 = ~ time +
    age +
    moves +
    cumulative_moves +
    cumulative_cost_lag +
    additional_cost +
    work_order_cost +
    box2 +
    age:I(cumulative_cost_lag + additional_cost) +
    age:work_order_cost,
  form3 = ~ time +
    ns(age, df = 5) +
    moves +
    cumulative_moves +
    cumulative_cost_lag +
    additional_cost +
    work_order_cost +
    box2 +
    age:I(cumulative_cost_lag + additional_cost) +
    age:work_order_cost,
  form5 = ~ unit2 * .,
  form6 = ~ unit2 * (time +
                       age +
                       moves +
                       cumulative_moves +
                       cumulative_cost_lag +
                       additional_cost +
                       work_order_cost +
                       box2 +
                       age:I(cumulative_cost_lag + additional_cost) +
                       age:work_order_cost),
  form7 = ~ unit2 * (time +
                       ns(age, df = 5) +
                       moves +
                       cumulative_moves +
                       cumulative_cost_lag +
                       additional_cost +
                       work_order_cost +
                       box2 +
                       age:I(cumulative_cost_lag + additional_cost) +
                       age:work_order_cost)

)

blip_formulas_stage_k <- list(
  form4 = ~ time +
    ns(age, df = 5) +
    ns(moves, df = 5) +
    ns(cumulative_moves, df = 5) +
    ns(cumulative_cost_lag, df = 5) +
    ns(additional_cost, df = 5) +
    ns(work_order_cost, df = 5) +
    box2 +
    age:I(cumulative_cost_lag + additional_cost) +
    age:work_order_cost,
  form8 = ~ unit2 * (time +
                       ns(age, df = 5) +
                       ns(moves, df = 5) +
                       ns(cumulative_moves, df = 5) +
                       ns(cumulative_cost_lag, df = 5) +
                       ns(additional_cost, df = 5) +
                       ns(work_order_cost, df = 5) +
                       box2 +
                       age:I(cumulative_cost_lag + additional_cost) +
                       age:work_order_cost)
)
blip_formulas_stage_1 <- list(
  form4 = ~ time +
    ns(age, df = 5) +
    ns(moves, df = 5) +
    ns(cumulative_moves, df = 5) +
    ns(additional_cost, df = 5) +
    ns(work_order_cost, df = 5) +
    box2 +
    age:I(additional_cost) +
    age:work_order_cost,
  form8 = ~ unit2 * (time +
                       ns(age, df = 5) +
                       ns(moves, df = 5) +
                       ns(cumulative_moves, df = 5) +
                       ns(additional_cost, df = 5) +
                       ns(work_order_cost, df = 5) +
                       box2 +
                       age:I(additional_cost) +
                       age:work_order_cost)
)

sl_blip_glm_stage_1 <- create.Learner("SL.glm_form",
                                      tune = list(formula = c(blip_formulas, blip_formulas_stage_1)),
                                      name_prefix = "blip_1")

sl_blip_glm_stage_k <- create.Learner("SL.glm_form",
                                      tune = list(formula = c(blip_formulas, blip_formulas_stage_k)),
                                      name_prefix = "blip_k")

# future arguments -------------------------------------------------------------

future.globals <- c(
  xgboost_learners_q$names,
  sl_blip_glm_stage_1$names,
  sl_blip_glm_stage_k$names
)

# simulation --------------------------------------------------------------

for (alpha in alpha0){
  print(alpha)
  # policy learner:
  pl <- policy_learn(
    type = "blip",
    control = control_blip(
      blip_models = c(
        q_sl(SL.library = sl_blip_glm_stage_1$names, env = environment()),
        replicate(K_partial0-1, q_sl(SL.library = sl_blip_glm_stage_k$names, env = environment()))
      )
    ),
    alpha = alpha
  )

  # simulation
  handlers(global = TRUE)
  plan(list(tweak("multisession", workers = workers)))
  res <- lava::sim(onerun_policy_eval,
                   R = R0,
                   args = list(
                     sim = sim,
                     n = n0,
                     n_approx = n_approx0,
                     K = K0,
                     K_partial = K_partial0,
                     q_models = q_models,
                     g_models = g_models,
                     policy_learn = pl,
                     M = M0,
                     save_g_functions = TRUE,
                     save_q_functions = TRUE
                   ),
                   future.packages = c("emrsim"),
                   future.globals = future.globals
  )

  out <- list(
    res = res,
    n = n0,
    n_approx = n_approx0,
    K = K0,
    K_partial = K_partial0,
    M = M0,
    alpha = alpha,
    summary_warnings = summary(warnings()),
  )

  # setting file name and saving
  tt <- Sys.Date()
  tt <- sub(" ", "_", tt)
  file_name <- paste(name,"_(alpha_", alpha,")_", tt, sep = "")
  file_name <- paste(file_name, ".rds", sep = "")
  saveRDS(out, file = file_name)
}
