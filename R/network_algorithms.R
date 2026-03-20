
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

get_beem_static_results <- function(bootstrap_results) {
  
  results_split <- bootstrap_results |>
    group_by(fraction, iteration) |>
    group_split()
  
  metrics <- lapply(results_split, function(df) {
    
    # Convert back to wide matrix
    b_est <- df |>
      select(taxon_from, taxon_to, weight) |>
      pivot_wider(names_from = taxon_to, values_from = weight) |>
      column_to_rownames("taxon_from") |>
      as.matrix()
    
    # Compute PEP
    pep <- seqtime::getPep(b_est)
    
    # Compute edge number
    values <- c(b_est[upper.tri(b_est)], b_est[lower.tri(b_est)])
    edge_num <- length(which(values != 0))
    
    data.frame(
      fraction  = unique(df$fraction),
      iteration = unique(df$iteration),
      pep       = pep,
      edge_num  = edge_num
    )
  })
  
  bind_rows(metrics)
}


# LIMITS

run_limits_bootstrap <- function(abundance_data, n_bootstrap = 500, subsample_pcts = c(0.4, 0.5, 0.6, 0.7, 0.8), max_attempts = 50, bagging.iter = 100, output_file = NULL) {
  
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
          limits(abundance_sub, bagging.iter = bagging.iter)$Aest,
          error = function(e) {
            cat("\n  Warning: iteration", iter, "failed:", conditionMessage(e), "\n")
            return(NULL)
          }
        )
        
        if (!is.null(network)) {
          iter <- iter + 1
          cat("\n  Bootstrap", iter, "/", n_bootstrap, "completed (", pct * 100, "% subsampling)\n")
          
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
  
  # Combine all results
  final_results <- bind_rows(all_results) |>
    select(fraction, iteration, taxon_from, taxon_to, weight)
  
  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    write.csv(final_results, output_file, row.names = FALSE)
    cat("Results saved to", output_file, "\n")
  }
  
  return(final_results)
  
}
  

get_limits_results <- function(bootstrap_results) {
  
  results_split <- bootstrap_results |>
    group_by(fraction, iteration) |>
    group_split()
  
  metrics <- lapply(results_split, function(df) {
    
    # Convert back to wide matrix
    b_est <- df |>
      select(taxon_from, taxon_to, weight) |>
      pivot_wider(names_from = taxon_to, values_from = weight) |>
      column_to_rownames("taxon_from") |>
      as.matrix()

    
    # Compute PEP
    pep <- seqtime::getPep(b_est)
    
    # Compute edge number
    values <- c(b_est[upper.tri(b_est)], b_est[lower.tri(b_est)])
    edge_num <- length(which(values != 0))
    
    data.frame(
      fraction = unique(df$fraction),
      iteration = unique(df$iteration),
      pep = pep,
      edge_num = edge_num
    )
  })
  
  bind_rows(metrics)
}


# FlashWeave

run_flashweave <- function(abundance_data, sensitive = TRUE, heterogeneous = FALSE) {
  
  julia_eval('import FlashWeave: learn_network, save_network, load_data')
  
  data_file <- tempfile(fileext = ".csv")
  netw_file <- tempfile(fileext = ".gml")
  write.csv(abundance_data, data_file, row.names = TRUE)
  
  data_file_julia <- gsub("\\\\", "/", data_file)
  netw_file_julia <- gsub("\\\\", "/", netw_file)
  
  julia_assign("data_file", data_file_julia)
  julia_assign("netw_file", netw_file_julia)
  julia_assign("sensitive", sensitive)
  julia_assign("heterogeneous", heterogeneous)
  
  julia_eval('net_result = learn_network(data_file, sensitive=sensitive, heterogeneous=heterogeneous, transposed=true)')
  julia_eval('save_network(netw_file, net_result)')
  
  # Read network with weights
  G <- igraph::read_graph(netw_file_julia, format = "gml")
  
  # Extract adjacency matrix with weights
  inferred_network <- igraph::as_adjacency_matrix(G, attr = "weight", sparse = FALSE)
  
  if (file.exists(data_file)) file.remove(data_file)
  if (file.exists(netw_file_julia)) file.remove(netw_file_julia)
  
  return(inferred_network)
}

run_flashweave_bootstrap <- function(abundance_data, n_bootstrap = 500, subsample_pcts = c(0.4, 0.5, 0.6, 0.7, 0.8), max_attempts = 50, sensitive = TRUE, heterogeneous = FALSE, output_file = NULL) {
  
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
          run_flashweave(abundance_sub, sensitive = sensitive, heterogeneous = heterogeneous),
          error = function(e) {
            cat("\n  Warning: iteration", iter, "failed:", conditionMessage(e), "\n")
            return(NULL)
          }
        )
        
        if (!is.null(network)) {
          iter <- iter + 1
          cat("\n  Bootstrap", iter, "/", n_bootstrap, "completed (", pct * 100, "% subsampling)\n")
          
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
  
  # Combine all results
  final_results <- bind_rows(all_results) |>
    select(fraction, iteration, taxon_from, taxon_to, weight)
  
  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    write.csv(final_results, output_file, row.names = FALSE)
    cat("Results saved to", output_file, "\n")
  }
  
  return(final_results)
}

get_flashweave_results <- function(bootstrap_results) {
  
  results_split <- bootstrap_results |>
    group_by(fraction, iteration) |>
    group_split()
  
  metrics <- lapply(results_split, function(df) {
    
    b_est <- df |>
      select(taxon_from, taxon_to, weight) |>
      pivot_wider(names_from = taxon_to, values_from = weight) |>
      column_to_rownames("taxon_from") |>
      as.matrix()
    
    # Compute PEP
    pep <- seqtime::getPep(b_est)
    
    # Compute edge number
    values <- c(b_est[upper.tri(b_est)], b_est[lower.tri(b_est)])
    edge_num <- length(which(values != 0))
    
    data.frame(
      fraction = unique(df$fraction),
      iteration = unique(df$iteration),
      pep = pep,
      edge_num = edge_num
    )
  })
  
  bind_rows(metrics)
}
