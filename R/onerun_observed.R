#' @export
onerun_observed <- function(sim,
                            n,
                            K){

  sim_data <- sim(
    n = n,
    K = K
  )
  pd <- polle::policy_data(
    type = "long",
    data = sim_data$stage_data,
    baseline_data = sim_data$baseline_data,
    id = "id",
    stage = "stage",
    action = "action",
    utility = "utility"
  )
  value <- mean(get_utility(pd)$U)
  value_sd <- sd(get_utility(pd)$U)


  out <- c(value = value,
           value_sd = value_sd)
  return(out)
}
