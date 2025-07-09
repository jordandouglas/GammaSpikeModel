#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)
library(coda)
library(dplyr)
library(tidyr)
library(tracerer)

# ---- Argument parsing ----
option_list <- list(
  make_option(c("--nsims"), type = "integer", help = "Number of replicates", metavar = "number"),
  make_option(c("--burnin"), type = "double", help = "Burn-in fraction (e.g., 0.2)", metavar = "fraction"),
  make_option(c("--out"), type = "character", default = "posterior_summary.rds",
              help = "Output file for serialized results", metavar = "file")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---- Validate required options ----
if (is.null(opt$nsims) || is.null(opt$burnin)) {
  stop("Error: --nsims and --burnin must be specified. Use --help for usage.")
}

nsims <- opt$nsims
burnin <- opt$burnin
out_file <- opt$out

# ---- Load true values ----
truth <- tracerer::parse_beast_tracelog_file("../truth/truth_param.log")
truth_trees <- tracerer::parse_beast_trees("../truth/truth_annotated.trees")

if (length(truth_trees) != nrow(truth)) {
  stop("Mismatch between number of trees and rows in truth table.")
}

tree_heights <- sapply(truth_trees, function(tr) max(ape::node.depth.edgelength(tr)))
tree_lengths <- sapply(truth_trees, function(tr) sum(tr$edge.length))

truth$tree.height <- tree_heights
truth$tree.treeLength <- tree_lengths

# Parameters to extract
params <- c("tree.height", "tree.treeLength", "stubs.nstubs", "nSampleAncestors",
            "netDiv", "lambda", "turnover", "samplingProportion", "rhoSampling",
            "clockMean", "clockSD", "spikeMean", "spikeShape", "gammaShape", "kappa")

# ---- Summarization Function ----
summarize_param <- function(param_name, nsims, burnin) {
  param_results <- vector("list", nsims)

  message(sprintf("Summarizing: %s", param_name))

  for (i in 1:nsims) {
    message(sprintf("  [Param: %s] Processing replicate %d / %d", param_name, i, nsims))
    path <- paste0("../templates/rep", i, "/estimatedParameters.log")

    if (!file.exists(path)) {
      message(sprintf("    Skipping: file not found"))
      next
    }

    post <- tryCatch({
      tracerer::parse_beast_tracelog_file(path)
    }, error = function(e) {
      message(sprintf("    Skipping: failed to parse (%s)", e$message))
      return(NULL)
    })

    if (is.null(post) || !(param_name %in% colnames(post))) {
      message("    Skipping: parameter missing in log")
      next
    }

    post_trimmed <- tracerer::remove_burn_ins(post, burn_in_fraction = burnin)
    post_vals <- post_trimmed[[param_name]]

    if (length(post_vals) <= 1 || all(is.na(post_vals))) {
      message("    Skipping: only one or no valid post-burnin sample")
      next
    }

    hpd <- HPDinterval(as.mcmc(post_vals))
    mean_val <- mean(post_vals)
    true_val <- truth[[param_name]][i]

    param_results[[i]] <- data.frame(
      replicate = i,
      true = true_val,
      mean = mean_val,
      hpd_low = hpd[1],
      hpd_high = hpd[2],
      in_hpd = true_val >= hpd[1] & true_val <= hpd[2]
    )
  }

  result_df <- do.call(rbind, param_results[!sapply(param_results, is.null)])
  return(list(data = result_df, n_plotted = nrow(result_df)))

}

# ---- Plot Function ----
plot_param <- function(result_list, param_name) {
  df <- result_list$data
  n_plotted <- result_list$n_plotted

  coverage <- mean(df$in_hpd, na.rm = TRUE)
  pearson <- cor(df$true, df$mean, use = "complete.obs")

  ggplot(df, aes(x = true, y = mean)) +
    geom_point(aes(color = in_hpd), size = 2) +
    geom_errorbar(aes(ymin = hpd_low, ymax = hpd_high), width = 0) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("TRUE" = "skyblue", "FALSE" = "tomato")) +
    labs(
      x = "Truth", y = "Estimated",
      title = param_name,
      subtitle = sprintf("Coverage = %.0f%%, Pearson = %.2f, N plooted = %d", coverage * 100, pearson, n_plotted)
    ) +
    theme_minimal()
}

# ---- Run summary and plotting ----
all_results <- list()
all_plots <- list()

for (param in params) {
  result_list <- summarize_param(param, nsims, burnin)
  plt <- plot_param(result_list, param)
  all_results[[param]] <- result_list$data
  all_plots[[param]] <- plt
}

# ---- Save results ----
saveRDS(list(data = all_results, plots = all_plots), file = out_file)
