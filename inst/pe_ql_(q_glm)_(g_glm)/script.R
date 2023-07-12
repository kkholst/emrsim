name <- "pe_ql_(q_glm)_(g_glm)"

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
R0 <- 1e2
workers <- 20
alpha0 <- c(0.04, 0.05, 0.075, 0.1)

# nuisance models:
q_models <- q_glm()
g_models <- g_glm()

# future
handlers(global = TRUE)
plan(list(tweak("multisession", workers = workers)))

set.seed(1)
for (alpha in alpha0){
  print(alpha)
  # policy learner:
  pl <- policy_learn(type = "ql",
                     alpha = alpha)

  # simulation

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
                   )
  )

  summary_warnings <- summary(warnings())

  out <- list(
    n = n0,
    n_approx = n_approx0,
    K = K0,
    K_partial = K_partial0,
    res = res,
    M = M0,
    alpha = alpha,
    summary_warnings =summary_warnings
  )

  # setting file name and saving
  tt <- Sys.Date()
  tt <- sub(" ", "_", tt)
  file_name <- paste(name,"_(alpha_", alpha,")_", tt, sep = "")
  file_name <- paste(file_name, ".rds", sep = "")
  saveRDS(out, file = file_name)
}
