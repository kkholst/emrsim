#' @export
onerun_policy_eval <- function(sim,
                               n,
                               n_approx,
                               K,
                               K_partial,
                               type = "dr",
                               M = 1,
                               policy_learn,
                               q_models = NULL,
                               g_models = NULL,
                               save_g_functions = TRUE,
                               save_q_functions = TRUE,
                               return_policy_eval = FALSE,
                               add_policy_eval = FALSE,
                               deterministic_rewards = FALSE,
                               action_prob_range = NULL,
                               scale = 1){

  # simulate policy data:
  sim_data <- sim(
    n = n,
    K = K,
    deterministic_rewards = deterministic_rewards,
    action_prob_range = action_prob_range,
    scale = scale
  )
  stage_data <- sim_data$stage_data
  baseline_data <- sim_data$baseline_data

  pd <- policy_data(
    type = "long",
    data = stage_data,
    baseline_data = baseline_data,
    id = "id",
    stage = "stage",
    action = "action",
    utility = "utility"
  )
  # partial policy data:
  pd <- partial(pd, K=K_partial)

  # policy evaluation:
  pe <- policy_eval(
    type = type,
    policy_learn = policy_learn,
    policy_data = pd,
    q_models = q_models,
    g_models = g_models,
    M = M,
    save_g_functions = save_g_functions,
    save_q_functions = save_q_functions
  )

  if(return_policy_eval)
    return(pe)

  # getting the g- and Q-function (only available if M = 1)
  gf <- get_g_functions(pe)
  qf <- get_q_functions(pe)

  # setting model statistic output:

  g_no_errors <- NA
  if(inherits(gf$all_stages$g_model$model, "SuperLearner")){
    g_no_errors <- as.numeric(!any(c(gf$all_stages$g_model$model$errorsInLibrary,
                                     gf$all_stages$g_model$model$errorsInCVLibrary)))
  }

  q_no_errors <- NA
  if(inherits(qf$stage_1$q_model$model, "SuperLearner")){
    q_no_errors <- as.numeric(!any(c(
      unlist(lapply(qf, function(x) x$q_model$model$errorsInLibrary)),
      unlist(lapply(qf, function(x) x$q_model$model$errorsInCVLibrary))
    )))
  }

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
  pf <- lapply(1:K_partial, function(k) get_policy_functions(po, stage = k))
  rm(po)

  if( (n_approx %% 1e4) != 0)
    stop("n_approx must be dividable by 1e4.")

  U <- c()
  for (j in 1:(n_approx/1e4)){
    sim_data_policy <- sim(
      n = 1e4,
      K = K,
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

  out <-c(
    value = mean(U),
    value_sd = sd(U),
    value_estimate = pe$value_estimate,
    value_estimate_vcov = value_estimate_vcov,
    value_estimate_ipw = pe$value_estimate_ipw,
    value_estimate_or = pe$value_estimate_or,
    g_no_errors = g_no_errors,
    q_no_errors = q_no_errors
  )
  if (add_policy_eval == TRUE){
    out <- append(out, list(pe = pe, pd = pd))
  }

  return(out)
}
