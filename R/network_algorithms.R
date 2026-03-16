
# # BEEM-Static
# 
# est <- beem2param(res_diseased)
# pep_inferred <- seqtime::getPep(est$b.est)
# 
# # non-zero entries in the upper and lower triangle (skipping diagonal) 
# values <- est$b.est[upper.tri(est$b.est)]
# values <- c(values,est$b.est[lower.tri(est$b.est)])
# 
# edge_num <- (length(which(values != 0))) 
# 
# cat("Number of edges in the inferred network:", edge_num, "\n")
# cat("Proportion of positive edges in the inferred network:", pep_inferred, "\n")
# 
# # Separate running into different cores in the computer (parallel computing) and set.seed() for reproducibility.
# 
# 
# library(parallel)
# library(doParallel)
# library(foreach)
# 
# numCores <- detectCores()
# cl <- makeCluster(numCores - 1)
# registerDoParallel(cl)
# 
# subsample_fracs <- c(0.6, 0.7, 0.8, 0.9)
# n_rounds <- 10
# 
# seeds <- expand_grid(frac = subsample_fracs, round = 1:n_rounds) |>
#   mutate(seed = row_number())
# 
# results <- foreach(
#   frac = rep(subsample_fracs, each = n_rounds),
#   round = rep(1:n_rounds, times = length(subsample_fracs)),
#   .packages = c("beemStatic"),
#   .combine = "rbind"
# ) %dopar% {
#   
#   healthy_sub <- healthy_df[, sample(ncol(healthy_df), floor(ncol(healthy_df) * frac))]
#   diseased_sub <- diseased_df[, sample(ncol(diseased_df), floor(ncol(diseased_df) * frac))]
#   
#   res_h <- tryCatch(beemStatic::func.EM(healthy_sub), error = function(e) NULL)
#   res_d <- tryCatch(beemStatic::func.EM(diseased_sub), error = function(e) NULL)
#   
#   healthy_pep <- if (!is.null(res_h)) beemStatic::getPep(beemStatic::beem2param(res_h)$b.est) else NA_real_
#   diseased_pep <- if (!is.null(res_d)) beemStatic::getPep(beemStatic::beem2param(res_d)$b.est) else NA_real_
#   
#   data.frame(frac = frac, round = round, healthy_pep = healthy_pep, diseased_pep = diseased_pep)
#   
# }
# 
# stopCluster(cl)
# 
# 
# subsample_fracs <- c(0.6, 0.7, 0.8, 0.9)
# n_rounds <- 10
# 
# seeds <- expand_grid(frac = subsample_fracs, round = 1:n_rounds) |>
#   mutate(seed = row_number())
# 
# results <- list()
# for (i in seq_len(nrow(seeds))) {
#   message(glue::glue("Round {i}/{nrow(seeds)}"))
#   set.seed(seeds$seed[i])
#   
#   healthy_sub <- healthy_df[, sample(ncol(healthy_df), floor(ncol(healthy_df) * seeds$frac[i]))]
#   diseased_sub <- diseased_df[, sample(ncol(diseased_df), floor(ncol(diseased_df) * seeds$frac[i]))]
#   
#   res_h <- tryCatch(beemStatic::func.EM(healthy_sub), error = function(e) NULL)
#   res_d <- tryCatch(beemStatic::func.EM(diseased_sub), error = function(e) NULL)
#   
#   results[[i]] <- data.frame(
#     frac = seeds$frac[i],
#     round = seeds$round[i],
#     healthy_pep = if (!is.null(res_h)) seqtime::getPep(beemStatic::beem2param(res_h)$b.est) else NA_real_,
#     diseased_pep = if (!is.null(res_d)) seqtime::getPep(beemStatic::beem2param(res_d)$b.est) else NA_real_
#   )
# }
# 
# results_df <- bind_rows(results)
