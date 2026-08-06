#!/usr/bin/env Rscript

## TriOmicNet main workflow
## Usage:
##   Rscript code/TriOmicNet.R
##   Rscript code/TriOmicNet.R BRCA
##   Rscript code/TriOmicNet.R LUAD
##   Rscript code/TriOmicNet.R PRAD
##   Rscript code/TriOmicNet.R ALL

required_packages <- c("readr", "igraph", "irlba")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them with:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
}

## ------------------------------
## Repository and argument setup
## ------------------------------
get_script_path <- function() {
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", full_args, value = TRUE)

  if (length(file_arg) > 0) {
    return(normalizePath(
      sub("^--file=", "", file_arg[1]),
      winslash = "/",
      mustWork = TRUE
    ))
  }

  ## Interactive fallback
  candidate <- file.path(getwd(), "code", "TriOmicNet.R")
  if (file.exists(candidate)) {
    return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
  }

  if (basename(getwd()) == "code") {
    candidate <- file.path(getwd(), "TriOmicNet.R")
    if (file.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }

  stop("Cannot determine the location of code/TriOmicNet.R.")
}

script_path <- get_script_path()
code_dir <- dirname(script_path)
repo_dir <- normalizePath(file.path(code_dir, ".."), winslash = "/", mustWork = TRUE)
data_dir <- file.path(repo_dir, "data")

if (!dir.exists(data_dir)) {
  stop("Data directory not found: ", data_dir)
}

## construct_layer.R may create files relative to the working directory.
setwd(repo_dir)

args <- commandArgs(trailingOnly = TRUE)
target <- if (length(args) >= 1) toupper(args[1]) else "BRCA"
allowed_targets <- c("BRCA", "LUAD", "PRAD", "ALL")

if (!target %in% allowed_targets) {
  stop("Argument must be one of: BRCA, LUAD, PRAD, ALL")
}

## ------------------------------
## General utility functions
## ------------------------------
assert_file <- function(path, label = "Required file") {
  if (!file.exists(path)) {
    stop(label, " not found: ", path)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

find_existing_file <- function(label, candidates) {
  candidates <- unique(candidates)
  hits <- candidates[file.exists(candidates)]

  if (length(hits) == 0) {
    stop(
      label, " was not found. Checked:\n",
      paste(candidates, collapse = "\n")
    )
  }

  if (length(hits) > 1) {
    message(label, ": multiple files found; using ", hits[1])
  }

  normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

load_rdata_env <- function(path) {
  path <- assert_file(path, "R data file")
  env <- new.env(parent = emptyenv())
  object_names <- load(path, envir = env)

  if (length(object_names) == 0) {
    stop("No objects were loaded from: ", path)
  }

  list(path = path, env = env, names = object_names)
}

select_loaded_object <- function(loaded, preferred_names, validator, label) {
  available <- loaded$names

  ## Exact preferred-name matching
  exact <- preferred_names[preferred_names %in% available]
  if (length(exact) > 0) {
    for (nm in exact) {
      obj <- loaded$env[[nm]]
      if (validator(obj)) return(obj)
    }
  }

  ## Case-insensitive preferred-name matching
  lower_available <- tolower(available)
  for (preferred in preferred_names) {
    idx <- which(lower_available == tolower(preferred))
    if (length(idx) > 0) {
      obj <- loaded$env[[available[idx[1]]]]
      if (validator(obj)) return(obj)
    }
  }

  ## Fall back only when exactly one object has the expected type
  valid_names <- available[vapply(
    available,
    function(nm) validator(loaded$env[[nm]]),
    logical(1)
  )]

  if (length(valid_names) == 1) {
    message(
      label, ": using object '", valid_names,
      "' loaded from ", loaded$path
    )
    return(loaded$env[[valid_names]])
  }

  stop(
    "Unable to identify ", label, " in ", loaded$path,
    "\nExpected object name(s): ", paste(preferred_names, collapse = ", "),
    "\nObjects actually present: ", paste(available, collapse = ", ")
  )
}

is_table_like <- function(x) {
  is.matrix(x) || is.data.frame(x)
}

as_numeric_matrix <- function(x, label) {
  if (!is_table_like(x)) {
    stop(label, " must be a matrix or data.frame.")
  }

  original_rownames <- rownames(x)
  original_colnames <- colnames(x)
  x_df <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)

  converted <- lapply(x_df, function(column) {
    suppressWarnings(as.numeric(as.character(column)))
  })
  out <- as.matrix(as.data.frame(converted, check.names = FALSE))

  rownames(out) <- original_rownames
  colnames(out) <- original_colnames

  if (any(!is.finite(out))) {
    stop(label, " contains non-numeric, NA, NaN, or Inf values.")
  }

  out
}

read_processed_expression <- function(path, expected_rows, expected_cols) {
  readers <- list(
    function() utils::read.table(
      path,
      sep = ",",
      header = FALSE,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    function() utils::read.table(
      path,
      sep = ",",
      header = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  )

  for (reader in readers) {
    raw <- tryCatch(reader(), error = function(e) NULL)
    if (is.null(raw)) next

    variants <- list(raw)
    if (ncol(raw) > 1) variants[[length(variants) + 1]] <- raw[, -1, drop = FALSE]
    if (nrow(raw) > 1) variants[[length(variants) + 1]] <- raw[-1, , drop = FALSE]
    if (nrow(raw) > 1 && ncol(raw) > 1) {
      variants[[length(variants) + 1]] <- raw[-1, -1, drop = FALSE]
    }

    for (candidate in variants) {
      if (nrow(candidate) != expected_rows || ncol(candidate) != expected_cols) next

      candidate_num <- suppressWarnings(as.matrix(data.frame(
        lapply(candidate, function(z) as.numeric(as.character(z))),
        check.names = FALSE
      )))

      if (all(is.finite(candidate_num))) return(candidate_num)
    }
  }

  stop(
    "The processed expression file has an incompatible format or dimension: ", path,
    "\nExpected dimensions: ", expected_rows, " rows x ", expected_cols, " columns."
  )
}

read_two_column_score <- function(path, label) {
  tab <- readr::read_csv(path, show_col_types = FALSE)
  if (ncol(tab) < 2) stop(label, " must contain at least two columns: Gene and Score.")

  out <- data.frame(
    Gene = as.character(tab[[1]]),
    Score = suppressWarnings(as.numeric(tab[[2]])),
    stringsAsFactors = FALSE
  )

  if (anyNA(out$Gene) || anyNA(out$Score)) {
    stop(label, " contains missing or non-numeric values.")
  }

  out
}

prepare_adjacency <- function(network, allowed_genes, label) {
  if (inherits(network, "igraph")) {
    network <- as.matrix(igraph::as_adjacency_matrix(network, sparse = FALSE))
  }

  if (!is_table_like(network)) {
    stop(label, " must be an adjacency matrix, data.frame, or igraph object.")
  }

  original_rownames <- rownames(network)
  original_colnames <- colnames(network)
  network <- as_numeric_matrix(network, label)
  rownames(network) <- original_rownames
  colnames(network) <- original_colnames

  if (nrow(network) != ncol(network)) {
    stop(label, " adjacency matrix is not square.")
  }
  if (is.null(rownames(network)) || is.null(colnames(network))) {
    stop(label, " adjacency matrix must have gene names as row and column names.")
  }

  common <- Reduce(intersect, list(
    rownames(network),
    colnames(network),
    as.character(allowed_genes)
  ))

  if (length(common) < 2) {
    stop(label, " has fewer than two genes overlapping the candidate gene list.")
  }

  network <- network[common, common, drop = FALSE]
  network[network != 0] <- 1
  diag(network) <- 0

  keep <- rowSums(network) > 0 | colSums(network) > 0
  network <- network[keep, keep, drop = FALSE]

  if (nrow(network) < 2) {
    stop(label, " contains no usable interactions after filtering.")
  }

  network
}

run_rwr <- function(network, seed_genes, label) {
  network_seeds <- intersect(as.character(seed_genes), rownames(network))

  if (length(network_seeds) == 0) {
    stop("No seed genes are present in the ", label, " network.")
  }

  result <- random_walk_ranking(network, network_seeds)
  result <- as.data.frame(result, stringsAsFactors = FALSE)

  if (ncol(result) < 2) {
    stop("random_walk_ranking() returned fewer than two columns for ", label, ".")
  }

  data.frame(
    Gene = as.character(result[[1]]),
    Score = suppressWarnings(as.numeric(result[[2]])),
    stringsAsFactors = FALSE
  )
}

map_scores <- function(target_genes, score_table) {
  idx <- match(target_genes, score_table$Gene)
  out <- numeric(length(target_genes))
  keep <- !is.na(idx)
  out[keep] <- score_table$Score[idx[keep]]
  out[!is.finite(out)] <- 0
  out
}

compute_sum_tni <- function(edge_index, number_of_genes) {
  neighbors <- vector("list", number_of_genes)

  if (nrow(edge_index) > 0) {
    for (r in seq_len(nrow(edge_index))) {
      i <- edge_index[r, 1]
      j <- edge_index[r, 2]
      neighbors[[i]] <- c(neighbors[[i]], j)
      neighbors[[j]] <- c(neighbors[[j]], i)
    }
  }

  neighbors <- lapply(neighbors, unique)
  degrees <- lengths(neighbors)
  sum_tni <- numeric(number_of_genes)

  if (nrow(edge_index) == 0) return(sum_tni)

  message("Calculating TNI on ", nrow(edge_index), " undirected STRING edges...")

  for (r in seq_len(nrow(edge_index))) {
    i <- edge_index[r, 1]
    j <- edge_index[r, 2]

    if (degrees[i] <= 1 || degrees[j] <= 1) next

    common_neighbors <- length(intersect(neighbors[[i]], neighbors[[j]]))
    edge_tni <- common_neighbors / min(degrees[i], degrees[j])

    sum_tni[i] <- sum_tni[i] + edge_tni
    sum_tni[j] <- sum_tni[j] + edge_tni
  }

  sum_tni
}

zscore <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (any(!is.finite(x))) stop("A score vector contains NA, NaN, or Inf values.")

  x_sd <- stats::sd(x)
  if (!is.finite(x_sd) || x_sd == 0) return(rep(0, length(x)))

  (x - mean(x)) / x_sd
}

## ------------------------------
## Main analysis for one cancer
## ------------------------------
run_triomiconet <- function(cancer) {
  prefix <- tolower(cancer)
  cancer_dir <- file.path(data_dir, cancer)
  result_dir <- file.path(repo_dir, "results", cancer)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(cancer_dir)) {
    stop("Cancer-specific data directory not found: ", cancer_dir)
  }

  message("\n========================================")
  message("Running TriOmicNet for ", cancer)
  message("Repository: ", repo_dir)
  message("========================================")

  helper_files <- c(
    DE_Score = file.path(code_dir, "DE_Score.R"),
    construct_layer = file.path(code_dir, "construct_layer.R"),
    random_walk = file.path(code_dir, "random_walk.R")
  )
  invisible(lapply(names(helper_files), function(label) {
    assert_file(helper_files[[label]], paste0(label, " script"))
  }))

  source(helper_files[["DE_Score"]])
  source(helper_files[["construct_layer"]])
  source(helper_files[["random_walk"]])

  if (!exists("DE_Score", mode = "function")) stop("DE_Score() was not defined by DE_Score.R.")
  if (!exists("construct_layer", mode = "function")) stop("construct_layer() was not defined by construct_layer.R.")
  if (!exists("random_walk_ranking", mode = "function")) {
    stop("random_walk_ranking() was not defined by random_walk.R.")
  }

  ## File names shown in the repository
  gene_exp_file <- file.path(cancer_dir, paste0(prefix, "_gene_exp.RData"))
  mut_file <- file.path(cancer_dir, paste0(prefix, "_mut.RData"))
  cancer_si_file <- file.path(cancer_dir, paste0(prefix, "_cancer_si.rda"))
  normal_si_file <- file.path(cancer_dir, paste0(prefix, "_normal_si.rda"))

  shared_files <- list(
    network = file.path(data_dir, "network.rda"),
    protein_info = file.path(data_dir, "protein_info.rda"),
    mir_gene_network = file.path(data_dir, "mir_gene_network.rda"),
    HINT = file.path(data_dir, "HINT.RData"),
    CPDB = file.path(data_dir, "CPDB.RData"),
    MULTINET = file.path(data_dir, "MULTINET.RData")
  )

  assert_file(gene_exp_file, "Gene-expression file")
  assert_file(mut_file, "Mutation file")
  assert_file(cancer_si_file, "Cancer miRNA file")
  assert_file(normal_si_file, "Normal miRNA file")
  invisible(lapply(names(shared_files), function(label) {
    assert_file(shared_files[[label]], paste0(label, " data file"))
  }))

  #### Step 1. Processing expression data ####
  message("Step 1/5: processing expression data")

  gene_exp_loaded <- load_rdata_env(gene_exp_file)
  gene_exp_obj <- select_loaded_object(
    gene_exp_loaded,
    preferred_names = c(paste0(prefix, "_gene_exp"), "mRNA_EXP", "gene_exp"),
    validator = is_table_like,
    label = paste0(cancer, " gene-expression matrix")
  )

  gene_exp_df <- as.data.frame(gene_exp_obj, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(gene_exp_df) < 2) {
    stop("The gene-expression matrix must contain a gene column and at least one sample column.")
  }

  gene_ids <- as.character(gene_exp_df[[1]])
  if (anyNA(gene_ids) || any(gene_ids == "")) stop("The gene column contains missing gene identifiers.")
  if (anyDuplicated(gene_ids)) stop("The gene-expression matrix contains duplicated gene identifiers.")

  mrna_exp <- suppressWarnings(as.matrix(data.frame(
    lapply(gene_exp_df[-1], function(z) as.numeric(as.character(z))),
    check.names = FALSE
  )))
  rownames(mrna_exp) <- gene_ids
  colnames(mrna_exp) <- colnames(gene_exp_df)[-1]

  if (any(!is.finite(mrna_exp))) {
    stop("The gene-expression matrix contains non-numeric, NA, NaN, or Inf values.")
  }

  means <- rowMeans(mrna_exp)
  vars <- apply(mrna_exp, 1, stats::var)
  sds <- apply(mrna_exp, 1, stats::sd)
  FF <- vars / (1 + vars)
  ET <- means + 2.5 * sds * FF

  gene_exp1 <- data.frame(
    gene = gene_ids,
    V1 = means,
    V2 = vars,
    V3 = sds,
    V4 = FF,
    V5 = ET,
    check.names = FALSE
  )

  utils::write.table(
    gene_exp1,
    file.path(result_dir, paste0("gene_exp1_", prefix, ".txt")),
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.csv(
    mrna_exp,
    file.path(result_dir, paste0("mrna_exp_", prefix, ".csv"))
  )
  utils::write.table(
    means,
    file.path(result_dir, paste0("means_", prefix, ".txt")),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  remrna_names <- c(
    paste0("remrna_exp_", prefix, ".txt"),
    paste0("remrna_exp_", prefix, ".txt.gz"),
    paste0("remrna_exp_", cancer, ".txt"),
    paste0("remrna_exp_", cancer, ".txt.gz")
  )
  remrna_dirs <- c(cancer_dir, code_dir, data_dir, repo_dir)
  remrna_file <- find_existing_file(
    paste0(cancer, " processed expression file"),
    unlist(lapply(remrna_dirs, function(d) file.path(d, remrna_names)))
  )

  remrna_exp <- read_processed_expression(
    remrna_file,
    expected_rows = nrow(mrna_exp),
    expected_cols = ncol(mrna_exp)
  )
  rownames(remrna_exp) <- rownames(mrna_exp)
  colnames(remrna_exp) <- colnames(mrna_exp)

  normal_idx <- grep("11A|11B", colnames(remrna_exp), ignore.case = TRUE)
  if (length(normal_idx) == 0) {
    stop("No normal samples containing 11A or 11B were found in the expression column names.")
  }

  cancer_idx <- setdiff(seq_len(ncol(remrna_exp)), normal_idx)
  if (length(cancer_idx) == 0) stop("No cancer expression samples remain after normal-sample selection.")

  normal_s <- remrna_exp[, normal_idx, drop = FALSE]
  cancer_s <- remrna_exp[, cancer_idx, drop = FALSE]

  cancer_si_loaded <- load_rdata_env(cancer_si_file)
  cancer_si_obj <- select_loaded_object(
    cancer_si_loaded,
    preferred_names = c(paste0(prefix, "_cancer_si"), "cancer_si"),
    validator = is_table_like,
    label = paste0(cancer, " cancer miRNA matrix")
  )

  normal_si_loaded <- load_rdata_env(normal_si_file)
  normal_si_obj <- select_loaded_object(
    normal_si_loaded,
    preferred_names = c(paste0(prefix, "_normal_si"), "normal_si"),
    validator = is_table_like,
    label = paste0(cancer, " normal miRNA matrix")
  )

  DE_05 <- DE_Score(normal_s, cancer_s, 0.5)
  mirna_value <- DE_Score(normal_si_obj, cancer_si_obj, 0.5)

  #### Step 2. Multilayer network construction ####
  message("Step 2/5: constructing the multilayer network")

  mut_loaded <- load_rdata_env(mut_file)
  mut_obj <- select_loaded_object(
    mut_loaded,
    preferred_names = c(
      paste0(prefix, "_mut"),
      paste0(prefix, "_G_mut"),
      "G_mut",
      "mutation"
    ),
    validator = is_table_like,
    label = paste0(cancer, " mutation matrix")
  )

  network_loaded <- load_rdata_env(shared_files$network)
  network_obj <- select_loaded_object(
    network_loaded,
    preferred_names = c("network"),
    validator = function(x) TRUE,
    label = "shared biological network"
  )

  protein_loaded <- load_rdata_env(shared_files$protein_info)
  protein_info_obj <- select_loaded_object(
    protein_loaded,
    preferred_names = c("protein_info"),
    validator = function(x) TRUE,
    label = "protein information"
  )

  ## This is the explicit Step 2 execution requested by the reviewer.
  construct_layer(mut_obj, network_obj, protein_info_obj, DE_05)

  generated_dirs <- unique(c(repo_dir, data_dir, cancer_dir, code_dir, result_dir))

  layer2_names <- c(
    paste0("layer_2_", prefix, ".csv"),
    paste0("layer_2_", cancer, ".csv"),
    paste0("layer2_", prefix, ".csv"),
    paste0("第二层所有节点集合_", prefix, ".csv"),
    "layer_2.csv"
  )
  c1_names <- c(paste0("C1_", prefix, ".csv"), paste0("C1_", cancer, ".csv"), "C1.csv")
  c2_names <- c(paste0("C2_", prefix, ".csv"), paste0("C2_", cancer, ".csv"), "C2.csv")
  string_names <- c(
    paste0("STRING_adj_", prefix, ".csv"),
    paste0("STRING_adj_", cancer, ".csv"),
    paste0("STRING要构成邻接矩阵的形式_", prefix, ".csv"),
    "STRING_adj.csv"
  )

  layer2_file <- find_existing_file(
    "Layer-2 node file generated by construct_layer()",
    unlist(lapply(generated_dirs, function(d) file.path(d, layer2_names)))
  )
  c1_file <- find_existing_file(
    "C1 file generated by construct_layer()",
    unlist(lapply(generated_dirs, function(d) file.path(d, c1_names)))
  )
  c2_file <- find_existing_file(
    "C2 file generated by construct_layer()",
    unlist(lapply(generated_dirs, function(d) file.path(d, c2_names)))
  )
  string_file <- find_existing_file(
    "STRING edge file generated by construct_layer()",
    unlist(lapply(generated_dirs, function(d) file.path(d, string_names)))
  )

  #### Step 3. Control ability and regulatory potential ####
  message("Step 3/5: calculating control ability and regulatory potential scores")

  mir_loaded <- load_rdata_env(shared_files$mir_gene_network)
  mir_gene_network_obj <- select_loaded_object(
    mir_loaded,
    preferred_names = c("mir_gene_network"),
    validator = is_table_like,
    label = "miRNA-gene network"
  )

  mir_gene_network <- as.matrix(mir_gene_network_obj)
  if (nrow(mir_gene_network) < 2 || ncol(mir_gene_network) < 2) {
    stop("mir_gene_network must contain a header row and at least two columns.")
  }
  colnames(mir_gene_network) <- as.character(mir_gene_network[1, ])
  mir_gene_network <- mir_gene_network[-1, , drop = FALSE]

  layer2 <- readr::read_csv(layer2_file, show_col_types = FALSE)
  if (ncol(layer2) < 1) stop("The layer-2 node file is empty.")

  candidate_genes <- as.character(layer2[[1]])
  if (anyNA(candidate_genes) || any(candidate_genes == "")) {
    stop("The layer-2 node file contains missing gene identifiers.")
  }

  two_gene_list <- data.frame(
    Gene = candidate_genes,
    S_Dir = 0,
    S_Ind = 0,
    `control ability score` = 0,
    TNI = 0,
    NL_TNI = 0,
    MIE = 0,
    MTCS = 0,
    STRING = 0,
    HINT = 0,
    CPDB = 0,
    MULTINET = 0,
    `multi-network diffusion score` = 0,
    `regulatory potential score` = 0,
    TriOmicNetScore = 0,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  mirna_value <- as.matrix(mirna_value)
  if (ncol(mirna_value) < 2) stop("DE_Score() returned an invalid miRNA score table.")

  for (i in seq_len(nrow(mirna_value))) {
    target_rows <- which(mir_gene_network[, 2] == mirna_value[i, 1])
    if (length(target_rows) == 0) next

    target_genes <- mir_gene_network[target_rows, 1]
    row_index <- match(target_genes, two_gene_list$Gene)
    row_index <- row_index[!is.na(row_index)]
    if (length(row_index) == 0) next

    increment <- suppressWarnings(as.numeric(mirna_value[i, 2]))
    if (!is.finite(increment)) next
    two_gene_list[row_index, "control ability score"] <-
      two_gene_list[row_index, "control ability score"] + increment
  }

  C1 <- read_two_column_score(c1_file, "C1")
  C2 <- read_two_column_score(c2_file, "C2")
  two_gene_list$S_Dir <- map_scores(two_gene_list$Gene, C1)
  two_gene_list$S_Ind <- map_scores(two_gene_list$Gene, C2)
  two_gene_list[["regulatory potential score"]] <- pmax(
    two_gene_list$S_Dir,
    two_gene_list$S_Ind
  )

  #### Step 4. Multi-network diffusion score ####
  message("Step 4/5: calculating the multi-network diffusion score")

  string_edges <- readr::read_csv(string_file, show_col_types = FALSE)
  if (ncol(string_edges) < 2) stop("The STRING edge file must contain at least two columns.")

  string_edges <- data.frame(
    GeneA = as.character(string_edges[[1]]),
    GeneB = as.character(string_edges[[2]]),
    stringsAsFactors = FALSE
  )
  string_edges <- string_edges[
    !is.na(string_edges$GeneA) &
      !is.na(string_edges$GeneB) &
      string_edges$GeneA != "" &
      string_edges$GeneB != "" &
      string_edges$GeneA != string_edges$GeneB,
    ,
    drop = FALSE
  ]

  genes <- unique(c(string_edges$GeneA, string_edges$GeneB))
  gene_index <- setNames(seq_along(genes), genes)
  edge_index <- cbind(
    as.integer(gene_index[string_edges$GeneA]),
    as.integer(gene_index[string_edges$GeneB])
  )
  edge_index <- cbind(
    pmin(edge_index[, 1], edge_index[, 2]),
    pmax(edge_index[, 1], edge_index[, 2])
  )
  edge_index <- unique(edge_index)

  N <- length(genes)
  PPI_STRING <- matrix(0, nrow = N, ncol = N, dimnames = list(genes, genes))
  PPI_STRING[edge_index] <- 1
  PPI_STRING[cbind(edge_index[, 2], edge_index[, 1])] <- 1

  ## Intended TNI calculation. This fixes the original ECC/TNI variable-name error
  ## and avoids the hard-coded 7686 x 7686 x 7686 loop.
  SUMTNI <- compute_sum_tni(edge_index, N)
  names(SUMTNI) <- genes
  tni_idx <- match(two_gene_list$Gene, names(SUMTNI))
  keep_tni <- !is.na(tni_idx)
  two_gene_list$TNI[keep_tni] <- SUMTNI[tni_idx[keep_tni]]

  maxTNI <- max(two_gene_list$TNI)
  if (!is.finite(maxTNI) || maxTNI == 0) maxTNI <- 1
  two_gene_list$NL_TNI <- two_gene_list$TNI / maxTNI

  mut_mat <- as_numeric_matrix(mut_obj, paste0(cancer, " mutation matrix"))
  Mutation <- t(mut_mat)

  if (is.null(rownames(Mutation))) {
    stop("The transposed mutation matrix must have gene names as row names.")
  }

  mutation_frequency <- rowSums(Mutation) / ncol(Mutation)
  vertex_W <- matrix(
    0.1,
    nrow = nrow(Mutation),
    ncol = ncol(Mutation),
    dimnames = dimnames(Mutation)
  )

  for (i in seq_len(ncol(Mutation))) {
    mutated_genes <- rownames(Mutation)[which(Mutation[, i] == 1)]
    valid_genes <- intersect(mutated_genes, two_gene_list$Gene)
    if (length(valid_genes) == 0) next

    p <- mutation_frequency[valid_genes]
    entropy_value <- ifelse(p > 0, -p * log2(p), 0)
    vertex_W[valid_genes, i] <- entropy_value
  }

  vertex_W[!is.finite(vertex_W)] <- 0
  MIE_all <- rowSums(vertex_W)
  mie_idx <- match(two_gene_list$Gene, names(MIE_all))
  keep_mie <- !is.na(mie_idx)
  two_gene_list$MIE[keep_mie] <- MIE_all[mie_idx[keep_mie]]
  two_gene_list$MTCS <- two_gene_list$NL_TNI * two_gene_list$MIE

  seed_table <- two_gene_list[order(-two_gene_list$MTCS), c("Gene", "MTCS")]
  top_n <- max(1, ceiling(0.05 * nrow(seed_table)))
  seed_genes <- seed_table$Gene[seq_len(top_n)]

  STRING_result <- run_rwr(PPI_STRING, seed_genes, "STRING")
  two_gene_list$STRING <- map_scores(two_gene_list$Gene, STRING_result)

  ## HINT can be stored either as mul_edge_list + mul_gene_list or as an adjacency matrix.
  hint_loaded <- load_rdata_env(shared_files$HINT)
  if (all(c("mul_edge_list", "mul_gene_list") %in% hint_loaded$names)) {
    mul_edge_list <- hint_loaded$env[["mul_edge_list"]]
    mul_gene_list <- hint_loaded$env[["mul_gene_list"]]

    from_idx <- suppressWarnings(as.integer(mul_edge_list[, 1]))
    to_idx <- suppressWarnings(as.integer(mul_edge_list[, 2]))
    valid_edge <-
      is.finite(from_idx) & is.finite(to_idx) &
      from_idx >= 1 & from_idx <= nrow(mul_gene_list) &
      to_idx >= 1 & to_idx <= nrow(mul_gene_list)

    hint_gene_edges <- cbind(
      as.character(mul_gene_list[from_idx[valid_edge], 2]),
      as.character(mul_gene_list[to_idx[valid_edge], 2])
    )
    hint_gene_edges <- hint_gene_edges[
      hint_gene_edges[, 1] != "" &
        hint_gene_edges[, 2] != "" &
        hint_gene_edges[, 1] != hint_gene_edges[, 2],
      ,
      drop = FALSE
    ]

    hint_graph <- igraph::graph_from_edgelist(hint_gene_edges, directed = FALSE)
    hint_genes <- intersect(igraph::V(hint_graph)$name, two_gene_list$Gene)
    if (length(hint_genes) < 2) stop("HINT has insufficient overlap with candidate genes.")
    hint_graph <- igraph::induced_subgraph(hint_graph, vids = hint_genes)
    HINT_mat <- as.matrix(igraph::as_adjacency_matrix(hint_graph, sparse = FALSE))
    HINT_mat <- prepare_adjacency(HINT_mat, two_gene_list$Gene, "HINT")
  } else {
    HINT_obj <- select_loaded_object(
      hint_loaded,
      preferred_names = c("HINT", "PPI"),
      validator = function(x) is_table_like(x) || inherits(x, "igraph"),
      label = "HINT network"
    )
    HINT_mat <- prepare_adjacency(HINT_obj, two_gene_list$Gene, "HINT")
  }

  HINT_result <- run_rwr(HINT_mat, seed_genes, "HINT")
  two_gene_list$HINT <- map_scores(two_gene_list$Gene, HINT_result)

  ## CPDB.RData and MULTINET.RData commonly both save the internal object as PPI.
  ## They are therefore loaded into separate environments to prevent overwriting.
  cpdb_loaded <- load_rdata_env(shared_files$CPDB)
  CPDB_obj <- select_loaded_object(
    cpdb_loaded,
    preferred_names = c("CPDB", "PPI"),
    validator = function(x) is_table_like(x) || inherits(x, "igraph"),
    label = "CPDB network"
  )
  CPDB_mat <- prepare_adjacency(CPDB_obj, two_gene_list$Gene, "CPDB")
  CPDB_result <- run_rwr(CPDB_mat, seed_genes, "CPDB")
  two_gene_list$CPDB <- map_scores(two_gene_list$Gene, CPDB_result)

  multinet_loaded <- load_rdata_env(shared_files$MULTINET)
  MULTINET_obj <- select_loaded_object(
    multinet_loaded,
    preferred_names = c("MULTINET", "PPI"),
    validator = function(x) is_table_like(x) || inherits(x, "igraph"),
    label = "MULTINET network"
  )
  MULTINET_mat <- prepare_adjacency(MULTINET_obj, two_gene_list$Gene, "MULTINET")
  MULTINET_result <- run_rwr(MULTINET_mat, seed_genes, "MULTINET")
  two_gene_list$MULTINET <- map_scores(two_gene_list$Gene, MULTINET_result)

  score_matrix <- as.matrix(two_gene_list[, c("STRING", "HINT", "CPDB", "MULTINET")])
  mean_score <- rowMeans(score_matrix)
  min_score <- apply(score_matrix, 1, min)
  max_score <- apply(score_matrix, 1, max)
  max_score[max_score == 0] <- 1e-6

  consistency_score <- min_score / max_score
  two_gene_list[["multi-network diffusion score"]] <- mean_score * consistency_score

  #### Step 5. SVD integration ####
  message("Step 5/5: integrating the three scores by SVD")

  score1_norm <- zscore(two_gene_list[["control ability score"]])
  score2_norm <- zscore(two_gene_list[["multi-network diffusion score"]])
  score3_norm <- zscore(two_gene_list[["regulatory potential score"]])
  S <- cbind(score1_norm, score2_norm, score3_norm)

  svd_result <- tryCatch(
    irlba::irlba(S, nv = 1, nu = 1),
    error = function(e) NULL
  )

  if (is.null(svd_result)) {
    message("irlba() could not be used; falling back to base R svd().")
    base_svd <- base::svd(S, nu = 1, nv = 1)
    U <- base_svd$u[, 1] * base_svd$d[1]
  } else {
    U <- svd_result$u[, 1] * svd_result$d[1]
  }

  two_gene_list$TriOmicNetScore <- abs(U)

  final_table <- two_gene_list[order(-two_gene_list$TriOmicNetScore), , drop = FALSE]
  final_table <- cbind(Rank = seq_len(nrow(final_table)), final_table)

  output_file <- file.path(result_dir, "final_score.csv")
  utils::write.csv(final_table, output_file, row.names = FALSE)

  message("Analysis finished for ", cancer)
  message("Final result: ", output_file)
  invisible(final_table)
}

cancers <- if (target == "ALL") c("BRCA", "LUAD", "PRAD") else target

for (cancer in cancers) {
  run_triomiconet(cancer)
  invisible(gc())
}
