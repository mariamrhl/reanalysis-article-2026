
plot_diversity_results <- function(df_res, cohort_name = "IBD") {
  
  df_res |>
    ggplot() +
    aes(x = group, y = pep, color = group, fill = group) +
    geom_violin(alpha = 0.4, color = NA) +
    geom_point(position = position_jitter(width = 0.1), alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Healthy" = "#0078B9", "Diseased" = "#EA0017")) +
    scale_fill_manual(values  = c("Healthy" = "#0078B9", "Diseased" = "#EA0017")) +
    facet_wrap(~ fraction) +
    labs(x = "", y = "Positive edge percentage (PEP)", color = "Group", fill = "Group",
         subtitle = paste("Comparison of PEP between healthy and diseased \nsamples in the", cohort_name, "cohort")) +
    stat_compare_means(
      comparisons = list(c("Healthy", "Diseased")),
      method = "wilcox.test",
      label = "p.signif",
      tip.length = 0.01,
      bracket.size = 0.8,
      fontface = "bold"
    ) +
    
    df_res |>
    ggplot() +
    aes(x = group, y = sheldon, color = group, fill = group) +
    geom_violin(alpha = 0.4, color = NA) +
    geom_point(position = position_jitter(width = 0.1), alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Healthy" = "#0078B9", "Diseased" = "#EA0017")) +
    scale_fill_manual(values  = c("Healthy" = "#0078B9", "Diseased" = "#EA0017")) +
    facet_wrap(~ fraction) +
    labs(x = "", y = "Sheldon evenness", color = "Group", fill = "Group",
         subtitle = paste("Comparison of evenness between healthy and \ndiseased samples in the", cohort_name, "cohort")) +
    stat_compare_means(
      comparisons = list(c("Healthy", "Diseased")),
      method = "wilcox.test",
      label = "p.signif",
      tip.length = 0.01,
      bracket.size = 0.8,
      fontface = "bold"
    ) +
    
    # shannon diversity
    df_res |>
    ggplot() +
    aes(x = group, y = shannon, color = group, fill = group) +
    geom_violin(alpha = 0.4, color = NA) +
    geom_point(position = position_jitter(width = 0.1), alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Healthy" = "#0078B9", "Diseased" = "#EA0017")) +
    scale_fill_manual(values  = c("Healthy" = "#0078B9", "Diseased" = "#EA0017")) +
    facet_wrap(~ fraction) +
    labs(x = "", y = "Shannon diversity", color = "Group", fill = "Group",
         subtitle = paste("Comparison of Shannon diversity between healthy and \ndiseased samples in the", cohort_name, "cohort")) +
    stat_compare_means(
      comparisons = list(c("Healthy", "Diseased")),
      method = "wilcox.test",
      label = "p.signif",
      tip.length = 0.01,
      bracket.size = 0.8,
      fontface = "bold"
    ) +
    
    # richness
    df_res |>
    ggplot() +
    aes(x = group, y = richness, color = group, fill = group) +
    geom_violin(alpha = 0.4, color = NA) +
    geom_point(position = position_jitter(width = 0.1), alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Healthy" = "#0078B9", "Diseased" = "#EA0017")) +
    scale_fill_manual(values  = c("Healthy" = "#0078B9", "Diseased" = "#EA0017")) +
    facet_wrap(~ fraction) +
    labs(x = "", y = "Richness", color = "Group", fill = "Group",
         subtitle = paste("Comparison of richness between healthy and \ndiseased samples in the", cohort_name, "cohort")) +
    stat_compare_means(
      comparisons = list(c("Healthy", "Diseased")),
      method = "wilcox.test",
      label = "p.signif",
      tip.length = 0.01,
      bracket.size = 0.8,
      fontface = "bold"
    ) +
    
    patchwork::plot_layout(nrow = 2)
}

plot_correlation_results <- function(df_res, variable_name) {
  
library(rlang)
df_res |>
  ggplot() +
  aes(x = sheldon, y = !!sym(variable_name), color = group, fill = group) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_smooth(aes(x = sheldon, y = !!sym(variable_name)), method = "lm",
              se = TRUE, inherit.aes = FALSE,
              color = "black", linetype = "dashed") +
  ggpubr::stat_cor(aes(label = after_stat(paste(r.label, p.label, sep = "~`,`~"))),
                   method = "pearson", label.x.npc = "left", label.y.npc = 0.12,
                   size = 3, show.legend = FALSE) +
  ggpubr::stat_cor(aes(x = sheldon, y = !!sym(variable_name)), method = "pearson",
                   inherit.aes = FALSE, color = "black",
                   label.x.npc = "left", label.y.npc = 0.015,
                   size = 3) +
  scale_color_manual(values = c("Healthy" = "#0078B9", "Diseased" = "#EA0017")) +
  scale_fill_manual(values  = c("Healthy" = "#0078B9", "Diseased" = "#EA0017")) +
  facet_wrap(~ fraction) 
  
}

plot_network_results <- function(df_res) {
  
  df_res |>
    ggplot() +
    aes(x = group, y = rho, color = group, fill = group) +
    geom_violin(alpha = 0.4, color = NA) +
    geom_point(position = position_jitter(width = 0.1), alpha = 0.6) +
    scale_color_manual(values = c("Healthy" = "#0078B9", "Diseased" = "#EA0017")) +
    scale_fill_manual(values  = c("Healthy" = "#0078B9", "Diseased" = "#EA0017")) +
    facet_wrap(~ fraction) +
    labs(x = "", y = "Net interactions", color = "Group", fill = "Group") +
    stat_compare_means(
      comparisons = list(c("Healthy", "Diseased")),
      method = "wilcox.test",
      label = "p.signif",
      tip.length = 0.01,
      bracket.size = 0.8,
      fontface = "bold"
    ) +
    
    df_res |>
    ggplot() +
    aes(x = group, y = pep, color = group, fill = group) +
    geom_violin(alpha = 0.4, color = NA) +
    geom_point(position = position_jitter(width = 0.1), alpha = 0.6) +
    scale_color_manual(values = c("Healthy" = "#0078B9", "Diseased" = "#EA0017")) +
    scale_fill_manual(values  = c("Healthy" = "#0078B9", "Diseased" = "#EA0017")) +
    facet_wrap(~ fraction) +
    labs(x = "", y = "Positive edge percentage (PEP)", color = "Group", fill = "Group") +
    stat_compare_means(
      comparisons = list(c("Healthy", "Diseased")),
      method = "wilcox.test",
      label = "p.signif",
      tip.length = 0.01,
      bracket.size = 0.8,
      fontface = "bold"
    ) +
    
    patchwork::plot_layout(nrow = 1)
  
}