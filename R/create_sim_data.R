#' @title simulate_sce: simulate a SingleCellExperiment as input data to COSICC_DA_group for testing and illustration
#' @description Generates simulated control and case `SingleCellExperiment` objects using the
#' `splatter` package. A binary column `tdTomato` is added to each cell to simulate a feature of interest, 
#' e.g. fluorescent marking of the cells.
#' 
#' @param nGenes Integer. Number of genes to simulate. Default is 10.
#' @param batchCells Integer. Number of cells to simulate.. Default is 5000.
#' @param nGroups Integer. Number of distinct cell types (groups). Default is 10.
#' @param seed Optional integer. Random seed for reproducibility. Default is `NULL`.
#'
#' @return A named list containing two `SingleCellExperiment` objects:
#' \describe{
#'   \item{case}{SCE object representing the case condition, including simulated depletion of tdTomato positive cells in some groups.}
#'   \item{control}{SCE object representing the control condition with randomly assigned tdTomato positive cells.}
#' }
#'
#' @importFrom splatter splatSimulate
#' @importFrom scater logNormCounts runPCA runUMAP
#' @importFrom SummarizedExperiment colData
#' @importFrom stats runif
#' @examples
#'   sim_data <- simulate_case_control_sce(seed = 123)
#'   sce_case <- sim_data$case
#'   sce_control <- sim_data$control
#' @export


simulate_sce <- function(
    nGenes = 10,
    batchCells = 5000,
    nGroups = 10,
    seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Simulate control data set
  sce_control <- splatter::splatSimulate(
    nGenes = nGenes,
    batchCells = batchCells,
    group.prob = rep(1 / nGroups, nGroups),
    method = "groups",
    verbose = FALSE
  )
  colnames(colData(sce_control))[colnames(colData(sce_control)) == "Group"] <- "celltype"
  
  sce_control$tdTomato <- "neg"
  for (j in seq_len(nGroups)) {
    group_idx <- which(sce_control$celltype == paste0("Group", j))
    sce_control$tdTomato[sample(group_idx, floor(length(group_idx) * 0.5 + 
                                                   (runif(1) < 0.5)))] <- "pos"
  }
  
  # Simulate case data set
  sce_case <- splatter::splatSimulate(
    nGenes = nGenes,
    batchCells = batchCells,
    group.prob = rep(1 / nGroups, nGroups),
    method = "groups",
    verbose = FALSE
  )
  colnames(colData(sce_case))[colnames(colData(sce_case)) == "Group"] <- "celltype"
  
  sce_case$tdTomato <- "neg"
  for (j in seq_len(nGroups)) {
    group_idx <- which(sce_case$celltype == paste0("Group", j))
    sce_case$tdTomato[sample(group_idx, floor(length(group_idx) * 0.5 + 
                                                (runif(1) <= 0.5)))] <- "pos"
  }
  
  # Deplete tdTomato+ cells from some groups
  tdTomato_cells_group1 <- which(sce_case$tdTomato == "pos" & sce_case$celltype == "Group1")
  tdTomato_cells_group2 <- which(sce_case$tdTomato == "pos" & sce_case$celltype == "Group2")
  tdTomato_cells_group3 <- which(sce_case$tdTomato == "pos" & sce_case$celltype == "Group3")
  
  depleted_cells <- c(
    sample(tdTomato_cells_group1, floor(0.5 * length(tdTomato_cells_group1) + (runif(1) < 0.5))),
    sample(tdTomato_cells_group2, floor(0.8 * length(tdTomato_cells_group2) + (runif(1) < 0.5))),
    sample(tdTomato_cells_group3, floor(0.95 * length(tdTomato_cells_group3) + (runif(1) < 0.5)))
  )
  
  sce_case <- sce_case[, -depleted_cells]
  
  # Simulating experimental bias
  marked_cells <- which(sce_case$tdTomato == "pos")
  unmarked_cells <- which(sce_case$tdTomato == "neg")
  sce_case <- 
    sce_case[,c(sample(marked_cells,length(marked_cells)),
                              sample(unmarked_cells,length(marked_cells)))]
  
  return(list(
    case = sce_case,
    control = sce_control
  ))
}



#' @title simulate_seurat: simulate a Seurat object as input data to 
#' COSICC_DA_group for testing and illustration
#'
#' @description Generates simulated control and case `Seurat` objects using the
#' `splatter` package. A binary column `tdTomato` is added to each cell to simulate a
#' feature of interest, e.g. fluorescent marking of the cells.
#'
#' @param nGenes Integer. Number of genes to simulate. Default is 10.
#' @param batchCells Integer. Number of cells to simulate. Default is 5000.
#' @param nGroups Integer. Number of distinct cell types (groups). Default is 10.
#' @param seed Optional integer. Random seed for reproducibility. Default is `NULL`.
#'
#' @return A named list containing two `Seurat` objects:
#' \describe{
#'   \item{case}{Seurat object representing the case condition, including simulated depletion of tdTomato-positive cells in some groups.}
#'   \item{control}{Seurat object representing the control condition with randomly assigned tdTomato-positive cells.}
#' }
#'
#' @importFrom splatter splatSimulate
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Seurat NormalizeData FindVariableFeatures ScaleData RunPCA RunUMAP
#' @importFrom stats runif
#' @examples
#'   sim_data <- simulate_seurat(seed = 123)
#'   seurat_case <- sim_data$case
#'   seurat_control <- sim_data$control
#' @export

simulate_seurat <- function(
    nGenes = 10,
    batchCells = 5000,
    nGroups = 10,
    seed = NULL
) {
  if (!requireNamespace("splatter", quietly = TRUE)) {
    stop("Package 'splatter' is required but not installed.")
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required but not installed.")
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  simulate_seurat_from_sce <- function(sce) {
    counts <- SummarizedExperiment::assay(sce, "counts")
    seurat_obj <- Seurat::CreateSeuratObject(counts = counts)
    seurat_obj$celltype <- sce$Group
    return(seurat_obj)
  }
  
  # Simulate control data
  sce_control <- splatter::splatSimulate(
    nGenes = nGenes,
    batchCells = batchCells,
    method = "groups",
    group.prob = rep(1 / nGroups, nGroups),
    verbose = FALSE
  )
  seurat_control <- simulate_seurat_from_sce(sce_control)
  seurat_control$tdTomato <- "neg"
  for (j in seq_len(nGroups)) {
    group_idx <- which(seurat_control$celltype == paste0("Group", j))
    seurat_control$tdTomato[sample(group_idx, floor(length(group_idx) * 0.5 +
                                                      (runif(1) < 0.5)))] <- "pos"
  }
  
  # Simulate case data
  sce_case <- splatter::splatSimulate(
    nGenes = nGenes,
    batchCells = batchCells,
    method = "groups",
    group.prob = rep(1 / nGroups, nGroups),
    verbose = FALSE
  )
  seurat_case <- simulate_seurat_from_sce(sce_case)
  seurat_case$tdTomato <- "neg"
  for (j in seq_len(nGroups)) {
    group_idx <- which(seurat_case$celltype == paste0("Group", j))
    seurat_case$tdTomato[sample(group_idx, floor(length(group_idx) * 0.5 +
                                                   (runif(1) <= 0.5)))] <- "pos"
  }
  
  # Simulate depletion in case
  group_filter <- function(seurat_obj, group_name) {
    which(seurat_obj$tdTomato == "pos" & seurat_obj$celltype == group_name)
  }
  
  tdTomato_cells_group1 <- group_filter(seurat_case, "Group1")
  tdTomato_cells_group2 <- group_filter(seurat_case, "Group2")
  tdTomato_cells_group3 <- group_filter(seurat_case, "Group3")
  
  depleted_cells <- c(
    sample(tdTomato_cells_group1, floor(0.5 * length(tdTomato_cells_group1) + (runif(1) < 0.5))),
    sample(tdTomato_cells_group2, floor(0.8 * length(tdTomato_cells_group2) + (runif(1) < 0.5))),
    sample(tdTomato_cells_group3, floor(0.95 * length(tdTomato_cells_group3) + (runif(1) < 0.5)))
  )
  
  seurat_case <- subset(seurat_case, cells = setdiff(colnames(seurat_case),
                                                     colnames(seurat_case)[depleted_cells]))
  marked_cells <- Seurat::WhichCells(seurat_case, expression = tdTomato == "pos")
  unmarked_cells <- Seurat::WhichCells(seurat_case, expression = tdTomato == "neg")
  unmarked_sub <- sample(unmarked_cells, length(marked_cells))
  
  # Subset Seurat object to balanced marked/unmarked cells
  seurat_case <- subset(seurat_case, cells = c(marked_cells, unmarked_sub))
  return(list(
    case = seurat_case,
    control = seurat_control
  ))
}

#' @title simulate_sce_kinetics: Simulate SingleCellExperiment objects with multiple 
#' lineaage trajectories and delay along some of the trajectories for perturbed cells
#' @description
#' Simulates control and case `SingleCellExperiment` objects using the 
#' splatter package's trajectory simulation.
#' Five lineage trajectories are simulated; for three of them, some tdTomato+ 
#' cells with higher pseudotime are depleted to simulate developmental delay, while for two, 
#' no such depletion is simulated.
#' The function allows for reproducibility via a user-specified random seed.
#'
#' @param batchCells Integer. Number of cells to simulate. Default is 500.
#' @param nGenes Integer. Number of genes to simulate. Default is 10.
#' @param de.prob Numeric. Probability of a gene being differentially expressed
#' along the lineage trajectory. 
#' Default is 0.5.
#' @param nPaths Integer. Number of trajectories (paths) to simulate. Default is 5.
#' @param seed Integer or NULL. Random seed for reproducibility. Default is NULL.
#' @param control_deplete_n Integer. Number of tdTomato+ high-pseudotime 
#' cells to deplete per lineage trajectory for the control data set. Default is 40.
#' @param case_deplete_n Integer. Number of tdTomato+ high-pseudotime cells 
#' to deplete per depleted trajectory for the case data set, e.g. the knockout chimera. Default is 80.
#' @param depleted_paths Character vector. Names of trajectories to deplete. 
#' Default is c("Path1", "Path2", "Path3").
#'
#' @return A named list containing two `SingleCellExperiment` objects:
#' \describe{
#'   \item{control}{SCE object representing the control condition. The \code{colData} contains:
#'     \itemize{
#'       \item \code{Cell}: Cell identifier
#'       \item \code{Path}: Lineage rajectory assignment (e.g., "Path1", "Path2", ...)
#'       \item \code{pt}: Pseudotime value for each cell
#'       \item \code{tdTomato}: Logical, whether the cell is tdTomato positive
#'       \item \code{nGene, nUMI, log10_total_counts, ...}: Additional QC columns from splatter
#'     }
#'   }
#'   \item{case}{SCE object representing the case condition, with additional depletion of tdTomato+ cells in selected trajectories. The \code{colData} columns are as above.}
#' }
#'
#' @importFrom splatter splatSimulatePaths
#' @importFrom scater logNormCounts
#' @import ggplot2
#' @export
simulate_sce_kinetics <- function(
    batchCells = 5000,
    nGenes = 10,
    de.prob = 0.5,
    nPaths = 5,
    seed = NULL,
    control_deplete_n = 10,
    case_deplete_n = "all", # can be "all" or a number
    depleted_paths = c("Path1", "Path2", "Path3")
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Simulate control with multiple trajectories
  sce_kinetics_control <- splatter::splatSimulatePaths(
    batchCells = batchCells,
    nGenes = nGenes,
    de.prob = de.prob,
    group.prob = rep(1/nPaths, nPaths),
    verbose = FALSE
  )
  
  # Rename Step to pt (pseudotime)
  sce_kinetics_control$pt <- sce_kinetics_control$Step
  sce_kinetics_control$Step <- NULL
  
  # Mark tdTomato+ cells (more marked in early pseudotime)
  n_tdTomato <- floor(0.5 * ncol(sce_kinetics_control))
  tdTomato_cells <- sample(seq_len(ncol(sce_kinetics_control)), n_tdTomato)
  sce_kinetics_control$tdTomato <- seq_len(ncol(sce_kinetics_control)) %in% tdTomato_cells
  
  # Deplete some tdTomato+ cells with higher pseudotime in selected lineage 
  # trajectories for the control data
  # (experimental rather than biological effect)
  to_keep <- rep(TRUE, ncol(sce_kinetics_control))
  for (path in depleted_paths) {
    idx_path <- which(sce_kinetics_control$Group == path)
    idx_tdTomato_highpt <- idx_path[sce_kinetics_control$tdTomato[idx_path] & 
                                      sce_kinetics_control$pt[idx_path] > 
                                      median(sce_kinetics_control$pt[idx_path])]
    n_deplete <- min(control_deplete_n, length(idx_tdTomato_highpt))
    if (n_deplete > 0) {
      to_deplete <- sample(idx_tdTomato_highpt, n_deplete)
      to_keep[to_deplete] <- FALSE
    }
  }
  sce_kinetics_control <- sce_kinetics_control[, to_keep]
  
  # Normalisation
  sce_kinetics_control <- scater::logNormCounts(sce_kinetics_control)
  
  # Simulate case (biological effect: remove only tdTomato+ high-pt cells in selected trajectories)
  to_keep_case <- rep(TRUE, ncol(sce_kinetics_control))
  for (path in depleted_paths) {
    idx_path <- which(sce_kinetics_control$Group == path)
    idx_tdTomato_highpt <- idx_path[sce_kinetics_control$tdTomato[idx_path] & 
                                      sce_kinetics_control$pt[idx_path] > median(sce_kinetics_control$pt[idx_path])]
    if (identical(case_deplete_n, "all")) {
      # Remove all tdTomato+ high-pt cells (only these)
      if (length(idx_tdTomato_highpt) > 0) {
        to_keep_case[idx_tdTomato_highpt] <- FALSE
      }
    } else {
      # Remove up to case_deplete_n tdTomato+ high-pt cells
      n_deplete <- min(case_deplete_n, length(idx_tdTomato_highpt))
      if (n_deplete > 0) {
        to_deplete <- sample(idx_tdTomato_highpt, n_deplete)
        to_keep_case[to_deplete] <- FALSE
      }
    }
  }
  sce_kinetics_case <- sce_kinetics_control[, to_keep_case]
return(list(
  control = sce_kinetics_control,
  case = sce_kinetics_case
))
}
