#globalVariables(c("cell_type", "sig"))
#' @import SingleCellExperiment
#' @import ggplot2
#' @import splatter
#' @import dplyr
#' @importFrom ggrepel geom_text_repel
#' @importFrom scales pseudo_log_trans
#' @importFrom stats median
#' @importFrom stats p.adjust
#' @importFrom stats fisher.test

#' @title COSICC_DA_group: Differential abundance analysis for groups of cells (e.g. cell types)
#' @description Performs differential abundance (DA) testing comparing marked versus unmarked cells 
#' between case and control single-cell data sets. For example, the case data set may be chimera data set with
#' a knock-out, and the control data set a wild-type chimera data set. This function supports both 
#' \code{SingleCellExperiment} and \code{Seurat} objects as inputs.
#'
#' The method subsamples cells to correct for experimental bias and uses 
#' Fisher's exact test to estimate the odds ratio and p-values 
#' for enrichment or depletion of cell types in marked relative to unmarked cells.
#' The function output includes a visualisation of the odds ratios and significance levels.
#'
#' @param sce_case A \code{SingleCellExperiment} or \code{Seurat} object for the case condition, e.g. a knock-out
#' chimera data set.
#' @param sce_control A \code{SingleCellExperiment} or \code{Seurat} object for the control condition, e.g. 
#' a wild-type chimera data set.
#' @param alpha Numeric scalar, significance cutoff for FDR (default 0.1).
#' @param thresh_marked_control Integer scalar, minimum number of marked cells in control to retain 
#' a cell type in the analysis (default 5).
#' @param thresh_unmarked_case Integer scalar, minimum number of unmarked cells in case to retain a cell 
#' type in the analysis (default 5).
#' @param marked_col Character, name of the meta.data (SeuratObject) or colData (SingleCellObject)
#' column indicating marked status (logical TRUE/FALSE) (default "marked").
#' @param celltype_col Character, name of the meta.data (SeuratObject) or colData (SingleCellObject)
#' column indicating cell types (default "cell_type").
#' @param seed Integer or NULL (default), allows the user to set a seed in the simulations for 
#' reproducibility 
#'
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{cell_type}{Cell type name.}
#'   \item{FDR}{FDR computed using the Benjamini-Hochberg method.}
#'   \item{odds_ratio}{Odds ratio of enrichment/depletion for each group of cells or cell type.}
#'   \item{sig}{Significance status: "enriched", "depleted", or "not significant".}
#' }
#' A plot showing the log10 odds ratio vs. -log10 FDR is also displayed. 
#'
#' @details
#' The function:
#' \itemize{
#'   \item Extracts and verifies metadata columns for marked status and cell types.
#'   \item Filters cell types based on minimum counts in marked control and unmarked case samples.
#'   \item Corrects for experimental bias in cell sampling by repeatedly subsampling to create synthetic
#' experiments with unbiased cell numbers.
#'   \item Produces a volcano-style plot highlighting enriched and depleted cell types.
#' }
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom stats fisher.test p.adjust median
#' @importFrom utils packageVersion
#' @import ggplot2 ggrepel scales
#'
#' @examples
#' \dontrun{
#' # Assuming sce_case and sce_control are Seurat or SingleCellExperiment objects with required metadata
#' result <- COSICC_DA_group(sce_case, sce_control)
#' head(result)
#' }
#'
#' @export

COSICC_DA_group <- function(
    sce_case, 
    sce_control, 
    alpha = 0.1, 
    thresh_marked_control = 5, 
    thresh_unmarked_case = 5,
    marked_col = "marked",
    celltype_col = "cell_type",
    seed = NULL
) {
  # Setting seed
  if (!(is.null(seed))){
    set.seed(seed)
  }
  
  # Helper function to retrieve necessary information from SingleCellExperiment or Seurat objects
  get_metadata <- function(obj, marked_col, celltype_col) {
    if (inherits(obj, "SingleCellExperiment")) {
      coldata <- SummarizedExperiment::colData(obj)
      md <- as.data.frame(coldata)
    } else if (inherits(obj, "Seurat")) {
      if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' required to handle Seurat objects, but it is not installed.")
      }
      md <- Seurat::FetchData(obj, vars = c(marked_col, celltype_col))
    } else {
      stop("Input object must be either a SingleCellExperiment or a Seurat object.")
    }
    
    # Check the input data
    if (!(marked_col %in% colnames(md))) {
      stop(paste0("Metadata must contain a '", marked_col, "' column with TRUE/FALSE values."))
    }
    if (!(celltype_col %in% colnames(md))) {
      stop(paste0("Metadata must contain a '", celltype_col, "' column."))
    }
    
    # Rename columns internally for ease of processing
    colnames(md)[colnames(md) == marked_col] <- "marked"
    colnames(md)[colnames(md) == celltype_col] <- "cell_type"
    
    # Check 'marked' is logical
    if (!is.logical(md$marked)) {
      stop("'marked' column must be logical (TRUE/FALSE).")
    }
    
    return(md)
  }
  
  # Extract information on cell type and marked status for each cell
  case_md <- get_metadata(sce_case, marked_col, celltype_col)
  control_md <- get_metadata(sce_control, marked_col, celltype_col)
  
  # Helper function to subset cells by marked status
  subset_cells <- function(obj, metadata, marked_flag) {
    cells <- rownames(metadata)[metadata$marked == marked_flag]
    if (inherits(obj, "SingleCellExperiment")) {
      return(obj[, cells])
    } else if (inherits(obj, "Seurat")) {
      return(Seurat::subset(obj, cells = cells))
    } else {
      stop("Input object must be either a SingleCellExperiment or a Seurat object.")
    }
  }
  
  sce_case_marked <- subset_cells(sce_case, case_md, TRUE)
  sce_case_unmarked <- subset_cells(sce_case, case_md, FALSE)
  sce_control_marked <- subset_cells(sce_control, control_md, TRUE)
  sce_control_unmarked <- subset_cells(sce_control, control_md, FALSE)
  
  # Extract cell_type vectors
  get_cell_types <- function(obj) {
    if (inherits(obj, "SingleCellExperiment")) {
      return(obj$cell_type)
    } else if (inherits(obj, "Seurat")) {
      return(Seurat::FetchData(obj, "cell_type")$cell_type)
    }
  }
  
  case_marked_ct <- get_cell_types(sce_case_marked)
  case_unmarked_ct <- get_cell_types(sce_case_unmarked)
  control_marked_ct <- get_cell_types(sce_control_marked)
  control_unmarked_ct <- get_cell_types(sce_control_unmarked)
  
  # Filter cell types by thresholds
  cell_types_control <- names(table(control_marked_ct))[table(control_marked_ct) >= thresh_marked_control]
  cell_types_case <- names(table(case_unmarked_ct))[table(case_unmarked_ct) >= thresh_unmarked_case]
  cell_types <- intersect(cell_types_control, cell_types_case)
  
  if (length(cell_types) == 0) {
    stop("No cell types pass the minimum cell count thresholds.")
  }
  
  # Compute ratios
  ratio_control <- sapply(cell_types, function(ct) {
    sum(control_marked_ct == ct) / sum(control_unmarked_ct == ct)
  })
  ratio_target <- sapply(cell_types, function(ct) {
    sum(case_marked_ct == ct) / sum(case_unmarked_ct == ct)
  })
  
  # Remove infinite or NA ratios
  valid_control <- ratio_control[!is.infinite(ratio_control) & !is.na(ratio_control)]
  valid_target <- ratio_target[!is.infinite(ratio_target) & !is.na(ratio_target)]
  
  norm_factor_control <- median(valid_control)
  norm_factor_target <- median(valid_target)
  
  # Helper function to obtain the overall number of cells
  n_cells <- function(obj) {
    if (inherits(obj, "SingleCellExperiment")) {
      return(ncol(obj))
    } else if (inherits(obj, "Seurat")) {
      return(ncol(obj@assays$RNA))
    }
  }
  
  sample_number_target <- 
    floor(min(n_cells(sce_case_marked), n_cells(sce_case_unmarked) / norm_factor_target))
  sample_number_control <- 
    floor(min(n_cells(sce_control_marked), n_cells(sce_control_unmarked) / norm_factor_control))
  
  # Subsampling and DA testing
  p_values <- vector("list", 100)
  odds_ratio <- vector("list", 100)
  
  for (k in seq_len(100)) {
    subsample_cells <- function(obj, n) {
      cells <- if (inherits(obj, "SingleCellExperiment")) colnames(obj) else colnames(obj@assays$RNA)
      sampled <- sample(cells, n)
      if (inherits(obj, "SingleCellExperiment")) {
        return(obj[, sampled])
      } else if (inherits(obj, "Seurat")) {
        return(Seurat::subset(obj, cells = sampled))
      }
    }
    
    case_marked_sub <- subsample_cells(sce_case_marked, sample_number_target)
    case_unmarked_sub <- 
      subsample_cells(sce_case_unmarked, floor(sample_number_target * norm_factor_target))
    control_marked_sub <- subsample_cells(sce_control_marked, sample_number_control)
    control_unmarked_sub <- 
      subsample_cells(sce_control_unmarked, floor(sample_number_control * norm_factor_control))
    
    
    case_marked_sub_ct <- if (inherits(case_marked_sub, "SingleCellExperiment")) 
      case_marked_sub$cell_type else Seurat::FetchData(case_marked_sub, "cell_type")$cell_type
    case_unmarked_sub_ct <- if (inherits(case_unmarked_sub, "SingleCellExperiment")) 
      case_unmarked_sub$cell_type else Seurat::FetchData(case_unmarked_sub, "cell_type")$cell_type
    control_marked_sub_ct <- if (inherits(control_marked_sub, "SingleCellExperiment")) 
      control_marked_sub$cell_type else Seurat::FetchData(control_marked_sub, "cell_type")$cell_type
    control_unmarked_sub_ct <- if (inherits(control_unmarked_sub, "SingleCellExperiment")) 
      control_unmarked_sub$cell_type else Seurat::FetchData(control_unmarked_sub, "cell_type")$cell_type
    
    odds_ratio[[k]] <- numeric(length(cell_types))
    p_values[[k]] <- numeric(length(cell_types))
    
    for (j in seq_along(cell_types)) {
      ct <- cell_types[j]
      contingency_table <- matrix(c(
        sum(case_marked_sub_ct == ct),
        sum(case_unmarked_sub_ct == ct),
        sum(control_marked_sub_ct == ct),
        sum(control_unmarked_sub_ct == ct)
      ), nrow = 2, byrow = TRUE)
      
      if (any(contingency_table == 0)) {
        odds_ratio[[k]][j] <- NA
        p_values[[k]][j] <- NA
      } else {
        test_res <- fisher.test(contingency_table)
        odds_ratio[[k]][j] <- test_res$estimate
        p_values[[k]][j] <- test_res$p.value
      }
    }
  }
  
  
  # Compute median p-values and odds ratios across the synthetic experiments 
  # created by subsampling
  p_values_mat <- do.call(cbind, p_values)
  odds_ratio_mat <- do.call(cbind, odds_ratio)
  median_p_values <- apply(p_values_mat, 1, median, na.rm = TRUE)
  median_odds_ratio <- apply(odds_ratio_mat, 1, median, na.rm = TRUE)
  
  FDR <- p.adjust(median_p_values, method = "BH")
  
  results_df <- data.frame(
    cell_type = cell_types,
    FDR = FDR,
    odds_ratio = median_odds_ratio,
    stringsAsFactors = FALSE
  )
  
  results_df <- results_df[order(results_df$FDR), ]
  results_df$sig <- "enriched"
  results_df$sig[results_df$odds_ratio < 1] <- "depleted"
  results_df$sig[results_df$FDR > alpha] <- "not significant"

  p <- ggplot(results_df, aes(x = log10(odds_ratio), y = -log10(FDR + 1e-100), text = cell_type)) +
    labs(
      x = expression(paste("log"[10], " of odds ratio")),
      y = expression(paste("-log"[10], " of FDR"))
    ) +
    geom_point(aes(color = sig), size = 3) +
    scale_color_manual(values = c(
      "not significant" = "grey",
      "enriched" = "darkblue",
      "depleted" = "darkred"
    ), name = "") +
    ggrepel::geom_text_repel(
      data = results_df[
        (apply(cbind(results_df$odds_ratio, 1 / results_df$odds_ratio), 1, max) > 1.5) |
          (results_df$p_values < 0.1),
      ],
      aes(label = cell_type),
      max.overlaps = Inf,
      size = 4
    ) +
    theme_classic(base_size = 14) +
    scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = 10)) +
    theme(
      axis.text = element_text(size = rel(0.75), color = "black"),
      axis.title = element_text(size = rel(1.0), color = "black")
    ) +
    annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = -Inf, fill = "blue", alpha = 0.1) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = Inf, ymax = -Inf, fill = "red", alpha = 0.1)
  
  print(p)
  
  return(results_df)
}

#' @title COSICC_DA_lineage
#' @details Computes differential cell fate probability for perturbation experiments
#' relative to control experiments.
#' @param sce_case SingleCellExperiment for the case data set,
#' colData needs to include the slots 'marked' (TRUE/FALSE) for the presence
#' of a perturbation (CRISPR, drug, disease etc.) or fluorescent marking,
#' and 'cell_type'
#' @param sce_control SingleCellExperiment for the control data set,
#' colData needs to include the slots 'marked' (TRUE/FALSE) for the presence
#' of a perturbation (CRISPR, drug, disease etc.)
#' @param lineage_scores data frame whose first column is the cell-id, and whose remaining columns
#' are scores reflecting the probability that a cell forms part of the lineage, e.g.
#' Waddington-OT cell fate matrices.
#' @param alpha significance level, default is 0.1
#' @return data.frame with the following columns: lineage, odds_ratio-odds ratio of the percentage of
#' marked positive cells among the target chimera cells belonging to the lineage versus
#' the respective ratio for the control experiment,
#' p_values, sig- whether a lineage is significantly enriched, depleted,
#' or neither significantly enriched nor depleted
#' @export

COSICC_DA_lineage <- function(sce_case,sce_control,lineage_scores,alpha=0.1){
  lineages <- colnames(lineage_scores)
  lineages <- setdiff(lineages,"id")
  lineage_scores$id <- rownames(lineage_scores)
  sce_case_marked <- sce_case[,as.logical(sce_case$marked)]
  sce_case_unmarked <- sce_case[,!(as.logical(sce_case$marked))]
  sce_control_unmarked <- sce_control[,!(as.logical(sce_control$marked))]
  sce_control_marked <- sce_control[,as.logical(sce_control$marked)]

  quantiles <- rep(NA,length(lineages))
  ratio_WT <- rep(NA,length(lineages))
  ratio_target <-rep(NA,length(lineages))

  names(ratio_target) <- lineages
  names(ratio_WT) <- lineages

  sce_case_score_marked <-
    lineage_scores[lineage_scores$id %in% sce_case_marked$cell,-ncol(lineage_scores)]
  sce_case_score_marked <-
    sce_case_score_marked[apply(sce_case_score_marked,1,function(x) max(x) > 10^-4),]


  sce_case_score_unmarked <-
    lineage_scores[lineage_scores$id %in% sce_case_unmarked$cell,colnames(lineage_scores)!="id"]
  sce_case_score_unmarked <-
    sce_case_score_unmarked[apply(sce_case_score_unmarked,1,function(x) max(x) > 10^-4),]


  sce_control_score <-
    lineage_scores[lineage_scores$id %in% sce_control_marked$cell,colnames(lineage_scores)!="id"]
  sce_control_score  <- sce_control_score[apply(sce_control_score,1,function(x) max(x) > 10^-4),]

  sce_control_score_neg <-
    lineage_scores[lineage_scores$id %in% sce_control_unmarked$cell,
                   colnames(lineage_scores)!="id"]
  sce_control_score_neg  <- sce_control_score_neg[apply(sce_control_score_neg ,1,
                                                  function(x) max(x) > 10^-4),]

  odds_ratio_all <- list()
  p_values_all <- list()

  for (jk in 1:30){
    sce_case_score_max <- apply(sce_case_score_marked,1,function(x) {y <- x>10^-4; temp <- x[y];
    return(sample(names(temp),1,prob=temp/sum(temp)))})
    #sce_case_score_max <- apply(sce_case_score_marked,1,function(x) return(names(x)[which.max(x)]))
    sce_case_score_max_table <- table(sce_case_score_max)
    sce_case_score_neg_max <- apply(sce_case_score_unmarked,1,function(x) {y <- x>10^-4; temp <- x[y];
    return(sample(names(temp),1,prob=temp/sum(temp)))})
    #sce_case_score_neg_max <- apply(sce_case_score_unmarked,1,function(x) return(names(x)[which.max(x)]))
    sce_case_score_neg_max_table <- table(sce_case_score_neg_max)
    sce_control_score_max <- apply(sce_control_score,1,function(x) {y <- x>10^-4; temp <- x[y];
    return(sample(names(temp),1,prob=temp/sum(temp)))})
    #sce_control_score_max <- apply(sce_control_score,1,function(x) return(names(x)[which.max(x)]))
    sce_control_score_max_table <- table(sce_control_score_max)
    sce_control_score_neg_max <- apply(sce_control_score_neg,1,function(x) {y <- x>10^-4; temp <- x[y];
    return(sample(names(temp),1,prob=temp/sum(temp)))})
    #sce_control_score_neg_max <- apply(sce_control_score_neg,1,function(x) return(names(x)[which.max(x)]))
    sce_control_score_neg_max_table <- table(sce_control_score_neg_max)


    nr_target_marked <- rep(0, length(lineages))
    names(nr_target_marked) <- lineages
    nr_target_marked[names(sce_case_score_max_table)] <- sce_case_score_max_table

    nr_target_unmarked <- rep(0, length(lineages))
    names(nr_target_unmarked) <- lineages
    nr_target_unmarked[names(sce_case_score_neg_max_table)] <- sce_case_score_neg_max_table

    nr_WT_marked <- rep(0, length(lineages))
    names(nr_WT_marked) <- lineages
    nr_WT_marked[names(sce_control_score_max_table)] <- sce_control_score_max_table

    nr_WT_unmarked <- rep(0, length(lineages))
    names(nr_WT_unmarked) <- lineages
    nr_WT_unmarked[names(sce_control_score_neg_max_table)] <- sce_control_score_neg_max_table

    xy <- nr_WT_marked >= 10
    lineages_keep <- names(nr_target_marked)[xy]

    nr_target_marked <- nr_target_marked[lineages_keep]
    nr_target_unmarked <- nr_target_unmarked[lineages_keep]
    nr_WT_marked <-  nr_WT_marked[lineages_keep]
    nr_WT_unmarked <- nr_WT_unmarked[lineages_keep]

    ratio_target <- nr_target_marked/nr_target_unmarked

    ratio_WT <- nr_WT_marked/nr_WT_unmarked

    norm_factor_target <- median(ratio_target)
    norm_factor_WT <- median(ratio_WT)

    sample_number_WT <- floor(min(length(sce_control_score_max ),length(sce_control_score_neg_max)/norm_factor_WT))
    sample_number_target <- floor(min(length(sce_case_score_max ),length(sce_case_score_neg_max)/norm_factor_target))
    odds_ratio <- list()
    p_values <- list()
    for (k in 1:30){
      sce_control_score_max <- sce_control_score_max[sample(1:length(sce_control_score_max),sample_number_WT)]
      sce_control_score_neg_max <- sce_control_score_neg_max[sample(1:length(sce_control_score_neg_max),sample_number_WT*norm_factor_WT)]

      sce_case_score_max <- sce_case_score_max[sample(1:length(sce_case_score_max),sample_number_target)]
      sce_case_score_neg_max <- sce_case_score_neg_max[sample(1:length(sce_case_score_neg_max),sample_number_target*norm_factor_target)]

      odds_ratio[[k]] <- rep(NA,length(lineages_keep))
      p_values[[k]] <- rep(NA,length(lineages_keep))

      for (j in 1:length(lineages_keep)){
        table_fisher <- c(sum(sce_case_score_max== lineages_keep[j]),sum(sce_case_score_neg_max == lineages_keep[j]),sum(sce_control_score_max == lineages_keep[j]),sum(sce_control_score_neg_max== lineages_keep[j]))
        dim(table_fisher) <- c(2,2)
        aa <- fisher.test(table_fisher)
        p_values[[k]][j] <- aa$p.value
        odds_ratio[[k]][j] <- aa$estimate
      }
    }
    odds_ratio <- do.call(cbind,odds_ratio)
    odds_ratio <- apply(odds_ratio,1,median)
    p_values <- do.call(cbind,p_values)
    p_values <- apply(p_values,1,median)
    p_values_all[[jk]] <- p_values
    odds_ratio_all[[jk]] <- odds_ratio
  }

  odds_ratio <- do.call(cbind,odds_ratio_all)
  odds_ratio <- apply(odds_ratio,1,median)
  p_values <- do.call(cbind,p_values_all)
  p_values <- apply(p_values,1,median)
  p_values <- p.adjust(p_values)

  fisher_test_lineages <- data.frame(lineage=lineages_keep,p_values=p_values,odds_ratio=odds_ratio)
  fisher_test_lineages <- fisher_test_lineages[order(fisher_test_lineages$p_values),]
  fisher_test_lineages$sig <- "enriched"
  fisher_test_lineages$sig[fisher_test_lineages$odds_ratio < 1] <- "depleted"
  fisher_test_lineages$sig[fisher_test_lineages$p_values> alpha] <- "not significant"


  p <- ggplot( fisher_test_lineages ,
    aes(x=log10(odds_ratio), y=-log10(p_values+1e-100),color=lineage)) +
    labs(x=expression(paste("log"[10]," of odds ratio")),
    y=expression(paste("-log"[10]," of FDR adjusted p-values"))) +
    geom_point(aes(color=sig),size=3) +
    scale_color_manual(values=c("not significant"="grey","enriched"="darkblue",
                                "depleted" = "darkred"),name="") +
    ggrepel::geom_text_repel(data=
             fisher_test_lineages[(apply(cbind(fisher_test_lineages$odds_ratio,1/fisher_test_lineages$odds_ratio),1,max)>1.5) |
                              (fisher_test_lineages$p_values < 0.1),],
                             aes(x=log10(odds_ratio), y=-log10(p_values+1e-100),
                                 label=lineage), max.overlaps=Inf, size=4) +
    theme_classic(base_size=14) +
    scale_y_continuous(trans=pseudo_log_trans(sigma=1,base=10))+
    theme(
      axis.text = element_text(size=rel(0.75), color='black'),
      axis.title = element_text(size=rel(1.0), color='black')) +
    annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = -Inf, fill= "blue",alpha=0.1) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = Inf, ymax = -Inf, fill= "red",alpha=0.1)
  print(p)
  return(fisher_test_lineages)
}
