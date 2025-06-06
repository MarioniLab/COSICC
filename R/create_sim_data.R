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