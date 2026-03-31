
# BEEM-Static

run_beem_static <- function(matrix, ncpu = 1, scaling = 1000, max_iter = 30, alpha = 1, lambda_choice = 1) {
  
  tryCatch({
    res <- beemStatic::func.EM(matrix, ncpu = ncpu, scaling = scaling, max.iter = max_iter, alpha = alpha, lambda.choice = lambda_choice)
    beemStatic::beem2param(res)$b.est}, 
    error = function(e) {cat("\nWarning: BEEM-Static failed:", conditionMessage(e), "\n")
      
    return(NULL)
  })
}


# LIMITS

run_limits <- function(matrix, bagging.iter = 100) {
  
  tryCatch({
    res <- limits(matrix, bagging.iter = bagging.iter)
    res$Aest}, 
    error = function(e) {cat("\nWarning: LIMITS failed:", conditionMessage(e), "\n")
      return(NULL)
    })
}


# FlashWeave

run_flashweave <- function(matrix, sensitive = TRUE, heterogeneous = FALSE) {
  
  data_file <- tempfile(fileext = ".csv")
  netw_file <- tempfile(fileext = ".gml")
  
  tryCatch({
    
    julia_eval('import FlashWeave: learn_network, save_network, load_data')
    
    write.csv(matrix, data_file, row.names = TRUE)
    
    data_file_julia <- gsub("\\\\", "/", data_file)
    netw_file_julia <- gsub("\\\\", "/", netw_file)
    
    julia_assign("data_file", data_file_julia)
    julia_assign("netw_file", netw_file_julia)
    julia_assign("sensitive", sensitive)
    julia_assign("heterogeneous", heterogeneous)
    
    julia_eval('net_result = learn_network(data_file, sensitive=sensitive, heterogeneous=heterogeneous, transposed=true)')
    julia_eval('save_network(netw_file, net_result)')
    
    G <- igraph::read_graph(netw_file_julia, format = "gml")
    igraph::as_adjacency_matrix(G, attr = "weight", sparse = FALSE)
    
  }, error = function(e) {
    cat("\nWarning: FlashWeave failed:", conditionMessage(e), "\n")
    return(NULL)
    
  }, finally = {
    if (file.exists(data_file)) file.remove(data_file)
    if (file.exists(netw_file)) file.remove(netw_file)
  })
}


run_network_bootstrap <- function(matrix, run_fn, n_bootstrap = 500, subsample_pcts = c(0.7, 0.8, 0.9), prevalence_threshold = 0.2, max_attempts = 50,
                          output_file = NULL, ...) {
  
  all_network_results   <- list()
  all_diversity_results <- list()
  
  for (pct in subsample_pcts) {
    cat("Running bootstrap with", pct * 100, "% subsampling\n")
    
    seen_hashes <- new.env(hash = TRUE, parent = emptyenv())
    iter     <- 0
    attempts <- 0
    
    while (iter < n_bootstrap && attempts < max_attempts) {
      attempts <- attempts + 1
      
      cat("\nBootstrap", iter, "/", n_bootstrap, "in progress (", pct * 100, "% subsampling)\n")
      
      matrix_sub <- subsample_unique_matrix(matrix, ncol(matrix), pct, seen_hashes)
      if (is.null(matrix_sub)) next
      
      matrix_filtered <- filter_prevalence_matrix(matrix_sub, prevalence_threshold)
      if (is.null(matrix_filtered) || nrow(matrix_filtered) < 5) next
      
      network <- run_fn(matrix_filtered, ...)
      if (is.null(network)) next
      
      iter <- iter + 1
      cat("\nBootstrap", iter, "/", n_bootstrap, "completed (", pct * 100, "% subsampling)\n")
      
      metrics <- compute_diversity_metrics(matrix_sub)
      
      all_network_results   <- c(all_network_results,   list(get_network_df(network, pct, iter, nrow(matrix_filtered))))
      all_diversity_results <- c(all_diversity_results, list(metrics |> mutate(fraction = pct, iteration = iter)))
    }
    
    if (attempts >= max_attempts) cat("\nWarning: max attempts reached at", pct * 100, "%\n")
    cat("\n", iter, "/", n_bootstrap, "successful iterations at", pct * 100, "%\n")
  }
  
  final_network <- dplyr::bind_rows(all_network_results)
  
  final_diversity_summary <- dplyr::bind_rows(all_diversity_results) |>
    group_by(fraction, iteration) |>
    summarise(across(c(shannon, richness, sheldon, pielou), mean), .groups = "drop")
  
  final_results <- final_network |>
    left_join(final_diversity_summary, by = join_by(fraction, iteration))
  
  if (!is.null(output_file)) {
    write.csv(final_results, output_file, row.names = FALSE)
    cat("Results saved to", output_file, "\n")
  }
  
  return(final_results)
}