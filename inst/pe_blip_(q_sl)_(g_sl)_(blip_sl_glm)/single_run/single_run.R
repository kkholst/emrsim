# Rscript -e "library(progressr)" -e "options(progressr.enable=TRUE)" -e "source('script.R')"

options(nwarnings = 1e4)

name <- "single_pe_alpha_(q_sl)_(g_sl)_(blip_sl_glm)"

library("polle")
library("emrsim")

set.seed(178)
n0 <- 2e4
K0 <- 16
n_approx0 <- 1e5
K_partial0 <- 8
M0 <- 1
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
  env = environment(),
  verbose = TRUE
)

# g-models ----------------------------------------------------------------

g_models <- g_sl(
  SL.library = c(
    "SL.mean",
    "SL.glm_form",
    "SL.gam_mgcv_g_2"
  ),
  cvControl = SuperLearner.CV.control(V = 2L),
  env = parent.frame(),
  verbose = TRUE
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

# simulate policy data:
sim_data <- sim(
  n = n0,
  K = K0
)

pd <- policy_data(
  type = "long",
  data = sim_data$stage_data,
  baseline_data = sim_data$baseline_data,
  id = "id",
  stage = "stage",
  action = "action",
  utility = "utility"
)
# partial policy data:
pd <- partial(pd, K=K_partial0)

res <- list()
for (a in seq_along(alpha0)){
  alpha <- alpha0[a]
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

  # policy evaluation:
  pe <- policy_eval(
    type = "dr",
    policy_learn = pl,
    policy_data = pd,
    q_models = q_models,
    g_models = g_models,
    M = M0
  )

  # getting the policy object from the policy evaluation if possible:
  po <- get_policy_object(pe)

  # fitting the policy on the full data
  if (is.null(po)){
    po <- policy_learn(
      policy_data = pd,
      g_models = g_models,
      q_models = q_models
    )
  }

  # simulating from the policy:
  pf <- lapply(1:K_partial0, function(k) get_policy_functions(po, stage = k))
  rm(po)

  if( (n_approx0 %% 1e4) != 0)
    stop("n_approx must be dividable by 1e4.")

  U <- c()
  for (j in 1:(n_approx0/1e4)){
    sim_data_policy <- sim(
      n = 1e4,
      K = K0,
      policy_functions = pf
    )
    pd_policy <- policy_data(
      type = "long",
      data = sim_data_policy$stage_data,
      baseline_data = sim_data_policy$baseline_data,
      id = "id",
      stage = "stage",
      action = "action",
      utility = "utility"
    )
    if (!all(get_action_set(pd_policy) == c("0", "1"))){
      stop("the action set is not 0/1")
    }
    U <- c(U, get_utility(pd_policy)$U)
  }

  # setting output:
  if (pe$type == "dr"){
    est <- estimate(pe)
    value_estimate_vcov <- vcov(est)
  } else{
    value_estimate_vcov <- NA
  }

  pe$g_functions <- NULL
  pe$q_functions <- NULL
  pe$policy_object <- NULL

  res[[a]] <- list(
    alpha = alpha,
    pe = pe,
    pd = pd,
    value = mean(U),
    value_sd = sd(U),
    value_estimate = pe$value_estimate,
    value_estimate_vcov = value_estimate_vcov,
    value_estimate_ipw = pe$value_estimate_ipw,
    value_estimate_or = pe$value_estimate_or
  )
}

# setting file name and saving
tt <- Sys.Date()
tt <- sub(" ", "_", tt)
file_name <- paste(name,"_", tt, sep = "")
file_name <- paste(file_name, ".rds", sep = "")
saveRDS(res, file = file_name)
