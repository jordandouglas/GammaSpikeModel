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
# truth_trees <- tracerer::parse_beast_trees("../truth/truth_annotated.trees")

# if (length(truth_trees) != nrow(truth)) {
#   stop("Mismatch between number of trees and rows in truth table.")
# }

# tree_heights <- sapply(truth_trees, function(tr) max(ape::node.depth.edgelength(tr)))
# tree_lengths <- sapply(truth_trees, function(tr) sum(tr$edge.length))

# truth$tree.height <- tree_heights
# truth$tree.treeLength <- tree_lengths

# ---- Parameters to extract ----
params <- c("tree.height", "tree.treeLength", "stubs.nstubs", "nSampledAncestors",
            "netDiv", "turnover", "samplingProportion", "lambda", "mu", "psi", "rhoSampling",
            "clockMean", "clockSD", "spikeMean", "spikeShape", "gammaShape", "kappa")

# ---- Efficient log reader ----
read_logs <- function(nsims) {
  log_data <- vector("list", nsims)

  for (i in 1:nsims) {
    log_data[[i]] <- list()
    log_main_path <- sprintf("../templates/rep%d/estimatedParameters.log", i)
    log_bdpsi_path <- sprintf("../templates/rep%d/estimatedParameters_bdpsi.log", i)

    log_data[[i]]$main <- if (file.exists(log_main_path)) {
      tryCatch(tracerer::parse_beast_tracelog_file(log_main_path), error = function(e) NULL)
    } else NULL

    log_data[[i]]$bdpsi <- if (file.exists(log_bdpsi_path)) {
      tryCatch(tracerer::parse_beast_tracelog_file(log_bdpsi_path), error = function(e) NULL)
    } else NULL
  }

  return(log_data)
}

# ---- Summarization Function ----
summarize_param <- function(param_name, log_data, truth, burnin) {
  nsims <- length(log_data)
  param_results <- vector("list", nsims)

  message(sprintf("Summarizing: %s", param_name))

  for (i in 1:nsims) {
    log_source <- if (param_name %in% c("lambda", "mu", "psi")) "bdpsi" else "main"
    post <- log_data[[i]][[log_source]]

    if (is.null(post)) {
      message(sprintf("  [rep %d] Skipping: %s log missing or unreadable", i, log_source))
      next
    }

    if (!(param_name %in% colnames(post))) {
      message(sprintf("  [rep %d] Skipping: parameter %s missing", i, param_name))
      next
    }

    post_trimmed <- tracerer::remove_burn_ins(post, burn_in_fraction = burnin)
    post_vals <- post_trimmed[[param_name]]

    if (length(post_vals) <= 1 || all(is.na(post_vals))) {
      message(sprintf("  [rep %d] Skipping: insufficient valid post-burnin samples", i))
      next
    }

    hpd <- HPDinterval(as.mcmc(post_vals))
    mean_val <- mean(post_vals)
    true_val <- truth[[param_name]][i]
    tolerance <- .Machine$double.eps^0.5
    in_hpd <- (true_val + tolerance) >= hpd[1] & (true_val - tolerance) <= hpd[2]

    param_results[[i]] <- data.frame(
      replicate = i,
      true = true_val,
      mean = mean_val,
      hpd_low = hpd[1],
      hpd_high = hpd[2],
      in_hpd = in_hpd
    #   in_hpd = true_val >= hpd[1] & true_val <= hpd[2]
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
      subtitle = sprintf("Coverage = %.0f%%, Pearson = %.2f, N plotted = %d", coverage * 100, pearson, n_plotted)
    ) +
    theme_minimal()
}

# ---- Run summaries and plots ----
log_data <- read_logs(nsims)

all_results <- list()
all_plots <- list()

for (param in params) {
  result_list <- summarize_param(param, log_data, truth, burnin)
  all_results[[param]] <- result_list$data
  all_plots[[param]] <- plot_param(result_list, param)
}

# ---- Save results ----
saveRDS(list(data = all_results, plots = all_plots), file = out_file)
