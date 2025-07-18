#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages({
  library(ape)
})

# Load parameter log
truth.df <- read.table("truth_param.log", header = TRUE, sep = "\t")
n <- nrow(truth.df)

# Load trees
all_trees <- read.nexus("truth.trees")
if (length(all_trees) != n) {
  stop("Mismatch: number of trees in 'truth_remaster.trees' != number of rows in 'truth_remaster_param.log'")
}

# Prepare output file
out_file <- file("truth_annotated.trees", open = "w")
cat("#NEXUS\n\nBegin trees;\n", file = out_file)

# Annotation function
annotate_tree_newick <- function(tree, tree_index, lambda, mu, psi, rho, clockRate, clockSD, spikeMean, spikeShape) {
  node_heights <- node.depth.edgelength(tree)
  tree_height <- max(node_heights)
  parents <- tree$edge[, 1]
  children <- tree$edge[, 2]
  edge_lengths <- tree$edge.length
  num_edges <- nrow(tree$edge)

  c_1 <- abs(sqrt((lambda - mu - psi)^2 + 4 * lambda * psi))
  c_2 <- -(lambda - mu - 2 * lambda * rho - psi) / c_1

  branches_df <- data.frame(
    edge_index = 1:num_edges,
    parent = parents,
    child = children,
    origin_MYA = tree_height - node_heights[parents],
    end_MYA = tree_height - node_heights[children],
    length = edge_lengths
  )

  # Add root edge if it exists
  has_root_edge <- !is.null(tree$root.edge)
  if (has_root_edge) {
    root_node <- Ntip(tree) + 1
    root_branch <- data.frame(
      edge_index = NA,
      parent = NA,
      child = root_node,
      origin_MYA = tree_height,
      end_MYA = tree_height - tree$root.edge,
      length = tree$root.edge
    )
    branches_df <- rbind(branches_df, root_branch)
  }

  delta_t <- with(branches_df, origin_MYA - end_MYA)
  numerator <- with(branches_df, (c_2 - 1) * exp(-c_1 * end_MYA) - (1 + c_2))
  denominator <- with(branches_df, (c_2 - 1) * exp(-c_1 * origin_MYA) - (1 + c_2))
  term1 <- delta_t * (lambda + mu + psi - c_1)
  term2 <- 2 * log(numerator / denominator)
  branches_df$expected_num_stubs <- term1 + term2

  if (has_root_edge) {
    branches_df$expected_num_stubs[nrow(branches_df)] <- 0
  }

  branches_df$num_stubs <- rpois(nrow(branches_df), lambda = branches_df$expected_num_stubs)
  if (has_root_edge) {
    branches_df$num_stubs[nrow(branches_df)] <- 0
  }

  branches_df$rel_clock_rate <- rlnorm(nrow(branches_df), meanlog = -(clockSD)^2 / 2, sdlog = clockSD)
  branches_df$abs_clock_rate <- branches_df$rel_clock_rate * clockRate
  branches_df$rel_spike_size <- rgamma(nrow(branches_df), shape = (branches_df$num_stubs + 1) * spikeShape, scale = 1 / spikeShape)
  branches_df$abs_spike_size <- branches_df$rel_spike_size * spikeMean

  build_newick <- function(node) {
    children <- which(parents == node)
    if (length(children) == 0) {
      label <- tree$tip.label[node]
    } else {
      subtrees <- sapply(tree$edge[which(tree$edge[, 1] == node), 2], build_newick)
      label <- paste0("(", paste(subtrees, collapse = ","), ")")
    }

    root_node <- Ntip(tree) + 1
    if (node == root_node && has_root_edge) {
      ann <- branches_df[nrow(branches_df), ]
      annotation <- sprintf("[&nstubs.branch=0.0,branchRates=%.16f,spikes=%.16f]",
                            ann$rel_clock_rate, ann$rel_spike_size)
      return(sprintf("%s%s:0.0", label, annotation))
    } else {
      edge_row <- which(tree$edge[, 2] == node)
      ann <- branches_df[edge_row, ]
      annotation <- sprintf("[&nstubs.branch=%.1f,branchRates=%.16f,spikes=%.16f]",
                            ann$num_stubs, ann$rel_clock_rate, ann$rel_spike_size)
      return(sprintf("%s%s:%.16f", label, annotation, ann$length))
    }
  }

  root_node <- Ntip(tree) + 1

  total_stubs <- sum(branches_df$num_stubs)
  return(list(
  newick = paste0(build_newick(root_node), ";"),
  total_stubs = total_stubs
  ))

}

# Initialize vector to store total stubs per tree
total_stubs_vec <- numeric(n)

# Annotate and write trees
for (i in 1:n) {
  tree <- all_trees[[i]]

  lambda <- truth.df[i, "lambda"]
  mu <- truth.df[i, "mu"]
  psi <- truth.df[i, "psi"]
  rho <- truth.df[i, "rhoSampling"]
  clockRate <- truth.df[i, "clockMean"]
  clockSD <- truth.df[i, "clockSD"]
  spikeMean <- truth.df[i, "spikeMean"]
  spikeShape <- truth.df[i, "spikeShape"]

  annotated_result <- annotate_tree_newick(tree, i - 1, lambda, mu, psi, rho, clockRate, clockSD, spikeMean, spikeShape)

  # Store total stubs
  total_stubs_vec[i] <- annotated_result$total_stubs

  #print(annotated)
  out_line <- paste0("tree STATE_", i - 1, " = ", annotated_result$newick)
  writeLines(out_line, con = out_file)
}

cat("End;\n", file = out_file)
close(out_file)

# Update and save truth_param.log
truth.df$stubs.nstubs <- total_stubs_vec
write.table(truth.df, file = "truth_param.log", sep = "\t", quote = FALSE, row.names = FALSE)

cat("Updated parameters written to truth_param.log\n")
cat("Annotated trees written to truth_annotated.trees\n")