name <- "benchmark"

library("polle")
library("emrsim")
library("future.apply")
library("progressr")

set.seed(1)
n0 <- 2e4
K0 <- 16
K_partial0 <- 8
R0 <- 1e2

handlers(global = TRUE)
plan("multisession")
res <- lava::sim(onerun_observed,
                 R = R0,
                 args = list(
                   sim = sim,
                   n = n0,
                   K = K0
                 )
)
plan("sequential")

# setting file name and saving
tt <- Sys.time()
tt <- sub(" ", "_", tt)
file_name <- paste(name, tt, sep = "_")
file_name <- paste(file_name, ".rds", sep = "")
saveRDS(res, file = file_name)





