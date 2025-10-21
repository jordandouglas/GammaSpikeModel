#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)
library(coda)
library(dplyr)
library(tidyr)
library(tracerer)

options(warn = 1) # Forces all warnings to print immediately

# ---- Argument parsing ----
option_list <- list(
  make_option(c("--nsims"), type = "integer", help = "Number of replicates", metavar = "number"),
  make_option(c("--burnin"), type = "double", help = "Burn-in fraction (e.g., 0.2)", metavar = "fraction"),
  make_option(c("--out"), type = "character", default = "posterior_summary",
              help = "Output file base name for results (.rds) and plots (.pdf)", metavar = "name")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---- Validate required options ----
if (is.null(opt$nsims) || is.null(opt$burnin)) {
  stop("Error: --nsims and --burnin must be specified. Use --help for usage.")
}

nsims <- opt$nsims
burnin <- opt$burnin
out_base_name <- opt$out

# ---- Load true values ----
truth_log_path <- "truth/truth_param.log"
if (!file.exists(truth_log_path)) {
  stop("Error: Main truth log file not found at: ", truth_log_path)
}
truth <- tracerer::parse_beast_tracelog_file(truth_log_path)

# ---- Automatically Discover Parameters ----
cols_to_exclude <- c("Sample", "posterior", "likelihood", "prior")
params_to_analyze <- setdiff(colnames(truth), cols_to_exclude)

message(sprintf("Found %d parameters to analyze from the truth log.", length(params_to_analyze)))
if (length(params_to_analyze) == 0) {
  stop("No valid parameters found in the truth log file.")
}


# ---- Function to process a SINGLE replicate (Memory-Efficient) ----
process_replicate <- function(i, truth_data, burnin_frac) {
  log_path <- sprintf("templates/rep%d/estimatedParameters.log", i)
  
  if (!file.exists(log_path)) {
    warning(sprintf("Skipping replicate %d: Log file not found.", i))
    return(NULL)
  }
  
  post <- tryCatch(tracerer::parse_beast_tracelog_file(log_path), error = function(e) NULL)
  
  if (is.null(post)) {
    warning(sprintf("Skipping replicate %d: Log file is unreadable.", i))
    return(NULL)
  }

  post_trimmed <- tracerer::remove_burn_ins(post, burn_in_fraction = burnin_frac)

  if (nrow(post_trimmed) <= 1) {
    warning(sprintf("Skipping replicate %d: insufficient samples.", i))
    return(NULL)
  }

  common_params <- intersect(params_to_analyze, colnames(post))
  
  results_list <- lapply(common_params, function(param_name) {
    post_vals <- post_trimmed[[param_name]]
    
    # Check for sufficient post-burnin values and print a warning if not.
    if (length(post_vals) <= 1 || all(is.na(post_vals))) {
      # warning(sprintf(
      #   "Skipping replicate %d for parameter '%s': insufficient post-burnin samples.",
      #   i,
      #   param_name
      # ))
      return(NULL)
    }

    mcmc_vals <- as.mcmc(post_vals)
    hpd <- HPDinterval(mcmc_vals)
    ess <- effectiveSize(mcmc_vals)
    
    mean_val <- mean(post_vals, na.rm = TRUE)
    true_val <- truth_data[[param_name]][i]
    
    tolerance <- .Machine$double.eps^0.5
    in_hpd <- (true_val + tolerance) >= hpd[1] & (true_val - tolerance) <= hpd[2]

    tibble(
      replicate = i,
      parameter = param_name,
      true = true_val,
      mean = mean_val,
      hpd_low = hpd[1],
      hpd_high = hpd[2],
      in_hpd = in_hpd,
      ess = ess
    )
  })
  
  bind_rows(results_list)
}

# ---- Main Processing Loop ----
all_results_list <- lapply(1:nsims, function(i) {
  message(sprintf("Processing replicate %d/%d...", i, nsims))
  process_replicate(i, truth_data = truth, burnin_frac = burnin)
})

all_results_df <- bind_rows(all_results_list)

if (nrow(all_results_df) == 0) {
  stop("Processing finished, but no valid data was found. Cannot generate plots.")
}

# ---- Save summary data to .rds file ----
rds_file_name <- paste0(out_base_name)
saveRDS(all_results_df, file = rds_file_name)
message(sprintf("\nSummary data saved to %s", rds_file_name))

# ---- Check for Poorly Estimated clockSD & spikeShape Params ----
message("\n--- Checking for Poorly Estimated Parameters ---")
params_to_check <- c("clockSD", "spikeShape")
failed_replicates_df <- all_results_df %>%
  filter(parameter %in% params_to_check & in_hpd == FALSE)

if (nrow(failed_replicates_df) > 0) {
  message("Found replicates where true value was not in the 95% HPD interval:")
  
  failed_replicates_df <- failed_replicates_df %>% arrange(parameter, replicate)
  
  for (i in 1:nrow(failed_replicates_df)) {
    row <- failed_replicates_df[i, ]
    cat(sprintf(
      "  - Rep: %-4d | Param: %-12s | True: %-8.4f | HPD: [%.4f, %.4f] | ESS: %.1f\n",
      row$replicate,
      row$parameter,
      row$true,
      row$hpd_low,
      row$hpd_high,
      row$ess
    ))
  }
} else {
  message("All replicates successfully captured the true values for 'clockSD' and 'spikeShape' within their 95% HPD intervals.")
}

# ---- Generate and save plots to a single PDF ----
pdf_file_name <- paste0(out_base_name, "_plots.pdf")
pdf(pdf_file_name, width = 8, height = 7)

message(sprintf("\nGenerating and saving plots to %s...", pdf_file_name))
all_params_found <- unique(all_results_df$parameter)

for (param in all_params_found) {
  df_param <- all_results_df %>% filter(parameter == param)
  
  if (nrow(df_param) == 0) next
  
  coverage <- mean(df_param$in_hpd, na.rm = TRUE)
  pearson <- cor(df_param$true, df_param$mean, use = "complete.obs")
  n_plotted <- nrow(df_param)

  p <- ggplot(df_param, aes(x = true, y = mean)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_errorbar(aes(ymin = hpd_low, ymax = hpd_high), width = 0, alpha = 0.3) +
    geom_point(aes(color = in_hpd), size = 2.5, alpha = 0.8) +
    scale_color_manual(values = c("TRUE" = "cornflowerblue", "FALSE" = "firebrick"), name = "True value\nin 95% HPD") +
    labs(
      x = "True Value", y = "Estimated Mean",
      title = param,
      subtitle = sprintf("Coverage = %.1f%%,  Pearson's r = %.3f,  N = %d", coverage * 100, pearson, n_plotted)
    ) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
  
  print(p)
  message(sprintf("  - Plotted '%s'", param))
}

dev.off()
message("Done.")