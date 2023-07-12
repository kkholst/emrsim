name <- "pe_ql_(q_sl)_(g_sl)"

library("polle")
library("emrsim")
library("future.apply")
library("progressr")

set.seed(1)
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

# future arguments -------------------------------------------------------------

future.globals <- c(
  xgboost_learners_q$names
)

# simulation --------------------------------------------------------------

set.seed(1)
for (alpha in alpha0){
  print(alpha)
  # policy learner:
  pl <- policy_learn(type = "ql",
                     alpha = alpha)

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
  print(summary(warnings()))

  out <- list(
    n = n0,
    n_approx = n_approx0,
    K = K0,
    K_partial = K_partial0,
    res = res,
    pl = pl,
    M = M0,
    alpha = alpha
  )

  # setting file name and saving
  tt <- Sys.Date()
  tt <- sub(" ", "_", tt)
  file_name <- paste(name,"_(alpha_", alpha,")_", tt, sep = "")
  file_name <- paste(file_name, ".rds", sep = "")
  saveRDS(out, file = file_name)
}
