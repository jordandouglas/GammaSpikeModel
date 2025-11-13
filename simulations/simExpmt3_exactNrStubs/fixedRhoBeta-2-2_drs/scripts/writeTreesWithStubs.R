#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages({
    library(ape)
    library(dplyr)
    library(phangorn)
    library(tidyr)
    library(treeio)
})

count_stubs <- function(tree, node, metadata, stub_counts, carry_stubs = 0) {
  children <- Descendants(tree, node, type = "children")
  
    # If node is a tip
  if (length(children) == 0) {
    tip_label <- tree$tip.label[node] # Get tip label
    edge_length <- tree$edge.length[which(tree$edge[,2] == node)] # Get edge length
    samp_type <- metadata$samp[match(tip_label, metadata$node)] # Get sampling type
    
    if (is.null(stub_counts[[tip_label]])) stub_counts[[tip_label]] <- 0 # Initialize if not present
    
    if (samp_type == "psiSample" && edge_length == 0) { # Zero-length psi sample
      stub_counts[[tip_label]] <- carry_stubs
    } else if (!is.na(samp_type) && samp_type != "death") {
      stub_counts[[tip_label]] <- stub_counts[[tip_label]] + carry_stubs
    }
    return(list(stub_counts = stub_counts))
  }
  
    # If node is internal
  left  <- children[1]
  right <- children[2]

  left_is_tip <- length(Descendants(tree, left, type = "children")) == 0
  right_is_tip <- length(Descendants(tree, right, type = "children")) == 0

  left_length <- tree$edge.length[which(tree$edge[,2] == left)]
  right_length <- tree$edge.length[which(tree$edge[,2] == right)]
  
  left_tips  <- tree$tip.label[Descendants(tree, left,  type = "tips")[[1]]]
  right_tips <- tree$tip.label[Descendants(tree, right, type = "tips")[[1]]]
  
  left_all_death  <- all(metadata$samp[match(left_tips,  metadata$node)] == "death", na.rm = TRUE)
  right_all_death <- all(metadata$samp[match(right_tips, metadata$node)] == "death", na.rm = TRUE)
  
  node_label <- ifelse(is.null(tree$node.label[node - Ntip(tree)]), 
                       as.character(node), 
                       tree$node.label[node - Ntip(tree)])
  if (is.null(stub_counts[[node_label]])) stub_counts[[node_label]] <- 0
  
    # Recursion
  if (left_all_death & !right_all_death) { # left all death
    if (right_is_tip && right_length == 0) {
      carry_stubs <- carry_stubs
      } else {
      carry_stubs <- carry_stubs + 1
      }
    res <- count_stubs(tree, right, metadata, stub_counts, carry_stubs)
    stub_counts <- res$stub_counts
    
  } else if (!left_all_death & right_all_death) { # right all death
    if (left_is_tip && left_length == 0) {
      carry_stubs <- carry_stubs
      } else {
      carry_stubs <- carry_stubs + 1
      }
    res <- count_stubs(tree, left, metadata, stub_counts, carry_stubs)
    stub_counts <- res$stub_counts
    
  } else if (!left_all_death & !right_all_death) { # both daughter lineages survive
    stub_counts[[node_label]] <- stub_counts[[node_label]] + carry_stubs
    carry_stubs <- 0
    res1 <- count_stubs(tree, left, metadata, stub_counts, carry_stubs)
    stub_counts <- res1$stub_counts
    res2 <- count_stubs(tree, right, metadata, stub_counts, carry_stubs)
    stub_counts <- res2$stub_counts
  } 
  # if both daughter lineages are dead, no recursion
  
  return(list(stub_counts = stub_counts))
}

count_stubs_tree <- function(tree, metadata) {
  stub_counts <- list()
  root_node <- Ntip(tree) + 1
  res <- count_stubs(tree, root_node, metadata, stub_counts, carry_stubs = 0)

  node_heights <- node.depth.edgelength(tree)
  tree_height <- max(node_heights)
  parents <- tree$edge[, 1]
  children <- tree$edge[, 2]
  edge_lengths <- tree$edge.length
  num_edges <- nrow(tree$edge)

  branch_df <- data.frame(
    edge_index = 1:num_edges,
    parent = parents,
    child = children,
    origin_MYA = tree_height - node_heights[parents],
    end_MYA = tree_height - node_heights[children],
    length = edge_lengths
  )

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
    branch_df <- rbind(branch_df, root_branch)
  }

  stub_tbl <- tibble(
    node = names(res$stub_counts),
    nstubs = as.numeric(res$stub_counts)
  )

  metadata_with_stubs <- metadata %>%
    left_join(stub_tbl, by = "node") %>%
    mutate(nstubs = replace_na(nstubs, 0)) %>%
    mutate(node = as.integer(node))
    
  metadata_with_stubs_join <- metadata_with_stubs %>%
    select(node, samp, nstubs)

  branch_df <- branch_df %>%
    left_join(metadata_with_stubs_join, by = c("child" = "node"))

  return(list(tree = tree, metadta = metadata_with_stubs, branch_df = branch_df))
}

collapse_binary_nodes <- function(df) {
  df$parent <- as.character(df$parent)
  df$child  <- as.character(df$child)

  repeat {
    parents <- unique(df$parent)
    parents <- parents[!is.na(parents)]   # skip root NA
    made_change <- FALSE

    for (parent in parents) {
      # only targets with exactly one downstream edge
      if (sum(df$parent == parent, na.rm = TRUE) != 1) next

      up_idx   <- which(df$child == parent)   # the edge that points into this node
      down_idx <- which(df$parent == parent)  # the (single) edge leaving this node

      # safety checks
      if (length(up_idx) != 1 || length(down_idx) != 1) {
        # unexpected topology: skip this parent for now
        next
      }

      # apply the transformation to the upstream edge (up_idx)
      df$child[up_idx]   <- df$child[down_idx]
      df$end_MYA[up_idx] <- df$end_MYA[down_idx]
      # new length = sum of the two edge lengths (you may prefer origin - end)
      df$length[up_idx]  <- df$length[up_idx] + df$length[down_idx]
      df$samp[up_idx]    <- df$samp[down_idx]
      df$nstubs[up_idx]  <- df$nstubs[down_idx]

      # remove the downstream row immediately
      df <- df[-down_idx, , drop = FALSE]

      # optional: reindex edge_index
      if ("edge_index" %in% names(df)) df$edge_index <- seq_len(nrow(df))

      made_change <- TRUE
      # break to recompute parents fresh because df changed
      break
    }

    if (!made_change) break
  }

  # reset row names
  rownames(df) <- NULL
  df
}

generate_pruned_newick <- function(branch_df, lambda, mu, psi, rho, clockRate, clockSD, spikeMean, spikeShape) {

  branch_df$parent <- as.integer(branch_df$parent)
  branch_df$child <- as.integer(branch_df$child)
  
    edges_df <- branch_df
  
  live_tip_ids <- edges_df %>% # Identify live tips (non-death samples)
      filter(!is.na(samp), samp != 'death') %>%
      pull(child) %>%
      unique()
  
  live_nodes <- unique(live_tip_ids) # Initialize live nodes with live tips
  parents_to_check <- live_tip_ids # Initialize with live tips
  
  while (length(parents_to_check) > 0) { # Traverse upwards to find all ancestors of live tips
      new_parents <- edges_df %>% # Find parents of current nodes to check
          filter(child %in% parents_to_check) %>%
          pull(parent) %>%
          unique()
      
      newly_live_parents <- setdiff(new_parents, live_nodes) # Identify newly found live parents
      
      live_nodes <- c(live_nodes, newly_live_parents) # Add newly found parents to live nodes
      parents_to_check <- newly_live_parents
  }
  
  pruned_edges_df <- edges_df %>% # Prune edges to include only live nodes
      filter(parent %in% live_nodes, child %in% live_nodes)

    pruned_edges_df <- collapse_binary_nodes(pruned_edges_df)

  c_1 <- abs(sqrt((lambda - mu - psi)^2 + 4 * lambda * psi))
  c_2 <- -(lambda - mu - 2 * lambda * rho - psi) / c_1

  delta_t <- with(pruned_edges_df, origin_MYA - end_MYA)
  numerator <- with(pruned_edges_df, (c_2 - 1) * exp(-c_1 * end_MYA) - (1 + c_2))
  denominator <- with(pruned_edges_df, (c_2 - 1) * exp(-c_1 * origin_MYA) - (1 + c_2))
  term1 <- delta_t * (lambda + mu + psi - c_1)
  term2 <- 2 * log(numerator / denominator)
  pruned_edges_df$expected_nstubs <- term1 + term2

  root_id <- pruned_edges_df %>%
      filter(is.na(parent)) %>%
      pull(child)
  root_row <- is.na(pruned_edges_df$parent)

    pruned_edges_df$length[root_row] <- 0 # Set root branch length to zero
    pruned_edges_df$nstubs[root_row] <- 0
    pruned_edges_df$expected_nstubs[root_row] <- 0
    # pruned_edges_df$rel_clock_rate[root_row] <- 0
    # pruned_edges_df$rel_spike_size[root_row] <- 0

  # Determine shape parameter for each branch based on parent's other branches
  shape_param <- numeric(nrow(pruned_edges_df))
  for (i in seq_len(nrow(pruned_edges_df))) {
    parent_i <- pruned_edges_df$parent[i]
    # Find all branches with the same parent (siblings, including self)
    sibling_rows <- which(pruned_edges_df$parent == parent_i)
    # Exclude self
    sibling_rows <- sibling_rows[sibling_rows != i]
    # Check if any sibling has length 0
    has_zero_length_sibling <- any(pruned_edges_df$length[sibling_rows] == 0)
    if (has_zero_length_sibling) {
      shape_param[i] <- (pruned_edges_df$nstubs[i] + 0) * spikeShape
    } else {
      shape_param[i] <- (pruned_edges_df$nstubs[i] + 1) * spikeShape
    }
  }

  pruned_edges_df$rel_clock_rate <- rlnorm(nrow(pruned_edges_df), meanlog = -(clockSD)^2 / 2, sdlog = clockSD)
  pruned_edges_df$rel_spike_size <- rgamma(nrow(pruned_edges_df), shape = shape_param, scale = 1 / spikeShape)

    tol <- 1e-12
    # Force zero for rows with length == 0 (sampled ancestors)
    zero_mask <- !is.na(pruned_edges_df$length) & abs(pruned_edges_df$length) <= tol & !is.na(pruned_edges_df$parent)
    # pruned_edges_df$rel_clock_rate[zero_mask] <- 0
    pruned_edges_df$rel_spike_size[zero_mask] <- 0

    # print(pruned_edges_df)

  if (nrow(pruned_edges_df) == 0 && length(live_nodes) <= 1) {
      return("();")
  }

    total_stubs <- sum(pruned_edges_df$nstubs[!is.na(pruned_edges_df$parent)], na.rm = TRUE) # Calculate total stubs, ignoring root edge
    total_expected_stubs <- sum(pruned_edges_df$expected_nstubs[!is.na(pruned_edges_df$parent)], na.rm = TRUE) # Calculate total expected stubs, ignoring root edge
  
  adj_list_pruned <- pruned_edges_df %>%
      group_by(parent) %>%
      summarise(children = list(child))
    # print(adj_list_pruned)
  adj_map_pruned <- setNames(adj_list_pruned$children, adj_list_pruned$parent)
  
  node_props_pruned <- pruned_edges_df %>% # Extract properties for pruned nodes
    select(child, length, samp, nstubs, expected_nstubs, rel_clock_rate, rel_spike_size) %>%
    setNames(c("id", "length", "samp", "nstubs", "expected_nstubs", "rel_clock_rate", "rel_spike_size"))
  
  props_map_pruned <- split(node_props_pruned, node_props_pruned$id) # Create a map from node ID to its properties

    df_to_newick_pruned <- function(node) {
    node_str <- as.character(node)
    
    props <- props_map_pruned[[node_str]]
    
    is_tip_pruned <- !(node_str %in% names(adj_map_pruned))
    
    if (is_tip_pruned && !is.na(props$samp)) {
      label <- paste0(props$samp, "_", node)
    } else {
      label <- "" 
    }
    
    # metadata <- sprintf("[&nstubs.branch=%d]", props$nstubs)
    metadata <- sprintf("[&nstubs.branch=%d,branchRates=%.16f,spikes=%.16f,expected.nstubs.branch=%.16f]",
                        props$nstubs, props$rel_clock_rate, props$rel_spike_size, props$expected_nstubs)
    
    branch_length <- sprintf(":%.16f", props$length)
    
    if (is_tip_pruned) {
      return(paste0(label, metadata, branch_length))
    } else {
      children <- adj_map_pruned[[node_str]]
      
      subtrees <- sapply(children, df_to_newick_pruned)
      
      children_str <- paste(subtrees, collapse = ",")
      
      return(paste0("(", children_str, ")", label, metadata, branch_length))
    }
    }

  newick_with_root_branch <- df_to_newick_pruned(root_id)
    # print(newick_with_root_branch)
  
  # last_paren_pos <- tail(gregexpr("\\)", newick_with_root_branch)[[1]], 1) # Find position of last closing parenthesis (metadata for root branch)
  
  # newick_string <- paste0(substr(newick_with_root_branch, 1, last_paren_pos), ";") # Remove metadata for root branch

  return(list(newick = paste0(newick_with_root_branch, ";"), total_stubs = total_stubs, total_expected_stubs = total_expected_stubs,  branch_df = pruned_edges_df))
}

# Load parameter log
truth_df <- read.table("truth_param.log", header = TRUE, sep = "\t")
n <- nrow(truth_df)

# Load trees
all_trees <- read.beast("truth-full.trees")
if (length(all_trees) != n) {
  stop("Mismatch: number of trees in 'truth_remaster.trees' != number of rows in 'truth_remaster_param.log'")
}

# Prepare output file
out_file <- file("truth-sampled-stubs.trees", open = "w")
cat("#NEXUS\n\nBegin trees;\n", file = out_file)

# Initialize vector to store total stubs per tree
total_stubs_vec <- numeric(n)
total_expected_stubs_vec <- numeric(n)
output_tree_height_vec <- numeric(n)
output_tree_length_vec <- numeric(n)

# Annotate and write trees
for (i in 1:n) {
  tree <- all_trees[[i]]

  lambda <- truth_df[i, "lambda"]
  mu <- truth_df[i, "mu"]
  psi <- truth_df[i, "psi"]
  rho <- truth_df[i, "rhoSampling"]
  clockRate <- truth_df[i, "clockMean"]
  clockSD <- truth_df[i, "clockSD"]
  spikeMean <- truth_df[i, "spikeMean"]
  spikeShape <- truth_df[i, "spikeShape"]

    result <- generate_pruned_newick(count_stubs_tree(tree@phylo, tree@data)$branch_df,
                                    lambda, mu, psi, rho, clockRate, clockSD, spikeMean, spikeShape)

  # Store total stubs
  total_stubs_vec[i] <- result$total_stubs
  total_expected_stubs_vec[i] <- result$total_expected_stubs

  output_tree <- tryCatch({
    read.tree(text = result$newick)
  }, error = function(e) {
    return(NULL)
  })
  if (is.null(output_tree)) {
    output_tree_height_vec[i] <- 0
    output_tree_length_vec[i] <- 0
  } else {
    output_tree_height_vec[i] <- max(node.depth.edgelength(output_tree))
    output_tree_length_vec[i] <- sum(output_tree$edge.length)
  }

  out_line <- paste0("tree STATE_", i - 1, " = ", result$newick)
  writeLines(out_line, con = out_file)
}

# Update and save truth_param.log
truth_df$stubs.nstubs <- total_stubs_vec
truth_df$stubs.expected.nstubs <- total_expected_stubs_vec
truth_df$tree.height <- output_tree_height_vec
truth_df$tree.treeLength <- output_tree_length_vec
write.table(truth_df, file = "truth_param.log", sep = "\t", quote = FALSE, row.names = FALSE)

cat("Parameters written to truth_param.log\n")
cat("Newick trees with stubs per branch written to truth-sampled-stubs.trees\n")