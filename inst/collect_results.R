library("data.table")

# Q-learning: glm ---------------------------------------------------------

files <- list.files("inst/pe_ql_(q_glm)_(g_glm)")
files <- files[grepl(".rds",files)]
files <- paste("inst/pe_ql_(q_glm)_(g_glm)/",files, sep = "")

tmp <- lapply(files, function(x) readRDS(x))
tab <- lapply(tmp, function(x) data.frame(
  alpha = x$alpha,
  value = x$res[,"value"],
  value_estimate = x$res[,"value_estimate"],
  value_estimate_vcov = x$res[, "value_estimate_vcov"],
  value_estimate_ipw = x$res[, "value_estimate_ipw"],
  value_estimate_or = x$res[, "value_estimate_or"]
))

tab <- as.data.table(do.call(what = "rbind", tab))

res <- cbind(policy_learner = "ql", tab, qm = "glm", gm = "glm", cf = NA)

# Q-learning: sl ---------------------------------------------------------

files <- list.files("inst/pe_ql_(q_sl)_(g_sl)")
files <- files[grepl(".rds",files)]
files <- paste("inst/pe_ql_(q_sl)_(g_sl)/",files, sep = "")

tmp <- lapply(files, function(x) readRDS(x))
tab <- lapply(tmp, function(x) data.frame(
  alpha = x$alpha,
  value = x$res[,"value"],
  value_estimate = x$res[,"value_estimate"],
  value_estimate_vcov = x$res[, "value_estimate_vcov"],
  value_estimate_ipw = x$res[, "value_estimate_ipw"],
  value_estimate_or = x$res[, "value_estimate_or"]
))

tab <- as.data.table(do.call(what = "rbind", tab))
tab <- cbind(policy_learner = "ql", tab, qm = "sl", gm = "gam", cf = NA)

res <- rbind(res, tab)

# DR blip Q mean ----------------------------------------------------------

files <- list.files("inst/pe_blip_(q_mean)_(g_sl)_(blip_sl_glm)")
files <- files[grepl(".rds",files)]
files <- paste("inst/pe_blip_(q_mean)_(g_sl)_(blip_sl_glm)/",files, sep = "")

tmp <- lapply(files, function(x) readRDS(x))
tab <- lapply(tmp, function(x) data.frame(
  alpha = x$alpha,
  value = x$res[,"value"],
  value_estimate = x$res[,"value_estimate"],
  value_estimate_vcov = x$res[, "value_estimate_vcov"],
  value_estimate_ipw = x$res[, "value_estimate_ipw"],
  value_estimate_or = x$res[, "value_estimate_or"]
))
tab <- as.data.table(do.call(what = "rbind", tab))
tab <- cbind(policy_learner = "blip", tab, qm = "mean", gm = "gam", cf = NA)

res <- rbind(res, tab, fill = TRUE)

# DR blip sl -----------------------------------------------------------------

files <- list.files("inst/pe_blip_(q_sl)_(g_sl)_(blip_sl_glm)")
files <- files[grepl(".rds",files)]
files <- paste("inst/pe_blip_(q_sl)_(g_sl)_(blip_sl_glm)/",files, sep = "")

tmp <- lapply(files, function(x) readRDS(x))
tab <- lapply(tmp, function(x) data.frame(
  alpha = x$alpha,
  # r = nrow(x$res),
  value = x$res[,"value"],
  value_estimate = x$res[,"value_estimate"],
  value_estimate_vcov = x$res[, "value_estimate_vcov"],
  value_estimate_ipw = x$res[, "value_estimate_ipw"],
  value_estimate_or = x$res[, "value_estimate_or"]
))
tab <- as.data.table(do.call(what = "rbind", tab))
tab <- cbind(policy_learner = "blip", tab, qm = "sl", gm = "gam", cf = NA)

res <- rbind(res, tab, fill = TRUE)


# saving ------------------------------------------------------------------

saveRDS(res, file = "inst/res_policy_evaluation.rds")

# view:

rm(tab)
tab <- readRDS("inst/res_policy_evaluation.rds")
setkey(tab, policy_learner, qm, alpha)
tab <- tab[!is.na(value)]

tab[, res_dr := value_estimate - value]
tab[, res_or := value_estimate_or - value]
tab[, res_ipw := value_estimate_ipw - value]

tab[, cov_lower := value_estimate - 1.96 * sqrt(value_estimate_vcov)]
tab[, cov_upper := value_estimate + 1.96 * sqrt(value_estimate_vcov)]

tab[, ind := (value <= cov_upper) & (value >= cov_lower)]

tab <- tab[, list(
  'r' = .N,
  'value' = mean(value),
  'SD (value)' = sd(value),
  'bias (DR)' = mean(res_dr),
  'RMSE (DR)' = sqrt(mean(res_dr^2)),
  'coverage (DR)' = mean(ind),
  'bias (IPW)' = mean(res_ipw),
  'RMSE (IPW)' = sqrt(mean(res_ipw^2)),
  'bias (OR)' = mean(res_or),
  'RMSE (OR)' = sqrt(mean(res_or^2))
),
by = list('pl' = policy_learner, qm, gm, alpha)
]

