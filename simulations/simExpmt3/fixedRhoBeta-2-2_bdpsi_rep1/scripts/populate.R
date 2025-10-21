if (!require(rjson , quietly = TRUE)) install.packages("rjson")
if (!require(ape , quietly = TRUE)) install.packages("ape")
library(rjson)
library(ape)

# ---- Parse Named Command-Line Arguments ----
args <- commandArgs(trailingOnly = TRUE)

# Default values
BURNIN <- 0
NSIMS <- 50

# Helper function to extract named args
get_arg_value <- function(name, default) {
  flag <- paste0("--", name)
  if (flag %in% args) {
    idx <- match(flag, args)
    return(args[idx + 1])
  } else {
    return(default)
  }
}

# Parse named arguments
BURNIN <- as.numeric(get_arg_value("burnin", BURNIN))
NSIMS <- as.integer(get_arg_value("nsims", NSIMS))

cat(sprintf("Using burnin = %s, nsims = %s\n", BURNIN, NSIMS))


# ---- Tree Extraction Function ----
getTrees = function(fileName) {
  file.in = readLines(fileName, warn = FALSE)
  trees = file.in[grep("tree STATE_", file.in)]
  trees = gsub(".+[=] ", "", trees)
  trees
}

# ---- Read Data ----
truth.df = read.table("truth/truth_param.log", header=TRUE, sep="\t")
species.trees = getTrees("truth/truth_annotated.trees")

# ---- Handle Burn-in and Sampling ----
burnin.start = floor(BURNIN * nrow(truth.df))
if (burnin.start <= 0) burnin.start = 1
include = seq(from = burnin.start, to = nrow(truth.df), length.out = NSIMS)
include = floor(include)

# ---- Process Simulations ----
for (simnum in 1:NSIMS) {
  rownum = include[simnum]
  f = paste0("templates/rep", simnum)
  cat(paste(f, "\n"))
  dir.create(f, showWarnings = FALSE)

  JSON = list()
  sub.df = truth.df[rownum,]

  for (x in colnames(sub.df)) {
    if (grepl("nstubs", x) || grepl("ntaxa", x) || grepl("nancestor", x)) next
    JSON[[x]] = sub.df[, x]
  }

  tree_string <- species.trees[rownum]
  JSON[["tree"]] = gsub("&", "&amp;", tree_string)
  write(JSON[["tree"]], paste0(f, "/trueTree.newick"))

  # Extract tip labels from Newick
  tree_obj <- read.tree(text = tree_string)
  JSON[["taxonRange"]] <- paste0(tree_obj$tip.label, collapse = ",")

  JSON_str = as.character(rjson::toJSON(JSON, indent = 1))
  write(JSON_str, paste0(f, "/var.json"))
}
