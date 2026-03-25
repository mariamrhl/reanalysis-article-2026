
# BEEM-Static function

run_beem_static <- function(abundance_data, ncpu = 1, scaling = 1000, max_iter = 30, alpha = 1, lambda_choice = 1) {
  
  res <- 
    beemStatic::func.EM(abundance_data, ncpu = ncpu, scaling = scaling, max.iter = max_iter, alpha = alpha, lambda.choice = lambda_choice)
  
  inferred_network <- beemStatic::beem2param(res)$b.est
  
  return(inferred_network)
}

run_beem_static_bootstrap <- function(abundance_data, n_bootstrap = 500, subsample_pcts = c(0.7, 0.8, 0.9), ncpu = 4, scaling = 1000, max_iter = 20, 
alpha = 1, lambda_choice = 1, max_attempts = 50, output_file = "results/beem_bootstrap.csv") {
  
  # Create output directory if it doesn't exist
  # dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  n_samples <- ncol(abundance_data)
  all_results <- list()
  
  for (pct in subsample_pcts) {
    cat("Running bootstrap with", pct * 100, "% subsampling\n")
    
    seen_hashes <- new.env(hash = TRUE, parent = emptyenv())
    iter <- 0
    attempts <- 0
    
    while (iter < n_bootstrap && attempts < max_attempts) {
      
      # Sample without replacement
      subsample_idx <- sample(seq_len(n_samples), 
                              size = floor(n_samples * pct), 
                              replace = FALSE)
      abundance_sub <- abundance_data[, subsample_idx]
      
      # Hash to avoid duplicate subsets
      subset_hash <- digest::digest(colnames(abundance_sub))
      
      if (!exists(subset_hash, envir = seen_hashes)) {
        assign(subset_hash, TRUE, envir = seen_hashes)
        
        network <- tryCatch(
          run_beem_static(abundance_sub, ncpu = ncpu, scaling = scaling, max_iter = max_iter, alpha = alpha, lambda_choice = lambda_choice),
          error = function(e) {
            cat("\n  Warning: iteration", iter, "failed:", conditionMessage(e), "\n")
            return(NULL)
          }
        )
        
        if (!is.null(network)) {
          iter <- iter + 1
          cat("\n  Bootstrap", iter, "/", n_bootstrap, "completed (", pct * 100, "% subsampling)\n")
          
          # Convert matrix to long format and append metadata
          network_long <- network |>
            as.data.frame() |>
            rownames_to_column("taxon_from") |>
            pivot_longer(-taxon_from, names_to = "taxon_to", values_to = "weight") |>
            mutate(fraction = pct, iteration = iter)
          
          all_results <- c(all_results, list(network_long))
        }
      }
      attempts <- attempts + 1
    }
    
    if (attempts >= max_attempts) {
      cat("\nWarning: max attempts reached at", pct * 100, "%\n")
    }
    
    cat("\n", iter, "/", n_bootstrap, "successful iterations at", pct * 100, "%\n")
  }
  
  # Combine all results and save to single CSV
  final_results <- bind_rows(all_results) |>
    select(fraction, iteration, taxon_from, taxon_to, weight)
  
  write.csv(final_results, output_file, row.names = FALSE)
  cat("Results saved to", output_file, "\n")
  
  return(final_results)
}

