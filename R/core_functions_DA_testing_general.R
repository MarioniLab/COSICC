globalVariables(c("cell_type", "sig"))
#' @import SingleCellExperiment
#' @import ggplot2
#' @import splatter
#' @import dplyr
#' @importFrom ggrepel geom_text_repel
#' @importFrom scales pseudo_log_trans
#' @importFrom stats median
#' @importFrom stats p.adjust
#' @importFrom stats fisher.test
#'
#' @title COSICC_DA_group
#' @details COSICC_DA_group performs differential abundance testing
#' for perturbation experiments relative to control experiments.
#' @param sce_case SingleCellExperiment for the case data set,
#' colData needs to include the slots 'marked' (TRUE/FALSE) for the presence
#' of a perturbation (CRISPR, drug, disease etc.) or fluorescent marking,
#' and 'cell_type'
#' @param sce_control SingleCellExperiment for the control data set,
#' colData needs to include the slots 'marked' (TRUE/FALSE) for the presence
#' of a perturbation (CRISPR, drug, disease etc.)
#' @param alpha significance level, default is 0.1
#' @param thresh_marked_control minimum number of marked cells per cell type in the control data set
#' @param thresh_unmarked_case minimum number of unmarked cells per cell type in the case data set
#' @return data.frame with the following columns: cell_type, odds_ratio-odds ratio of the percentage of
#' marked cells among the cells in the case data versus the respective ratio for the control experiment,
#' p_values, sig (whether a cell_type is significantly enriched, depleted, or not enriched)
#' @export
COSICC_DA_group <- function(sce_case,sce_control,alpha=0.1,thresh_marked_control=5,thresh_unmarked_case=5){
  sce_case_marked <- sce_case[,sce_case$marked]
  sce_case_unmarked <- sce_case[,!(sce_case$marked)]
  sce_control_marked <- sce_control[,sce_control$marked]
  sce_control_unmarked <- sce_control[,!(sce_control$marked)]
  cell_types <- names(table(sce_control_marked$cell_type))[table(sce_control_marked$cell_type) >= thresh_marked_control]
  cell_types_1 <- names(table(sce_case_unmarked$cell_type))[table(sce_case_unmarked$cell_type) >= thresh_unmarked_case]
  cell_types <- intersect(cell_types,cell_types_1)
  ratio_control <- rep(0,length(cell_types))
  ratio_target <- rep(0, length(cell_types))
  for (j in 1:length(cell_types)){
    ratio_control[j] <- sum(sce_control_marked$cell_type == cell_types[j])/sum(sce_control_unmarked$cell_type == cell_types[j])
    ratio_target[j] <- sum(sce_case_marked$cell_type == cell_types[j])/sum(sce_case_unmarked$cell_type == cell_types[j])
  }
  names(ratio_target) <- cell_types
  names(ratio_control) <- cell_types

  norm_factor_control <- median(ratio_control[(!(is.infinite(ratio_control)))&(!(is.na(ratio_control)))])
  norm_factor_target <- median(ratio_target[(!(is.infinite(ratio_target)))&(!(is.na(ratio_target)))])
  sample_number_target <- floor(min(ncol(sce_case_marked),ncol(sce_case_unmarked)/norm_factor_target))
  sample_number_control <- floor(min(ncol(sce_control_marked),ncol(sce_control_unmarked)/norm_factor_control))
  p_values <- list()
  odds_ratio <- list()
  for (k in 1:100){
    sce_case_marked <- sce_case_marked[,sample(1:ncol(sce_case_marked),sample_number_target)]
    sce_case_unmarked <- sce_case_unmarked[,sample(1:ncol(sce_case_unmarked),floor(sample_number_target*norm_factor_target))]

    sce_control_marked <- sce_control_marked[,sample(1:ncol(sce_control_marked),sample_number_control)]
    sce_control_unmarked <- sce_control_unmarked[,sample(1:ncol(sce_control_unmarked),floor(sample_number_control*norm_factor_control))]

    odds_ratio[[k]] <- rep(0,length(cell_types))
    p_values[[k]] <- rep(0,length(cell_types))
    for (j in 1:length(cell_types)){
      table_fisher <- c(sum(sce_case_marked$cell_type == cell_types[j]),sum(sce_case_unmarked$cell_type == cell_types[j]),
                        sum(sce_control_marked$cell_type == cell_types[j]),sum(sce_control_unmarked$cell_type == cell_types[j]))
      dim(table_fisher) <- c(2,2)
      aa <- fisher.test(table_fisher)
      odds_ratio[[k]][j] <- aa$estimate
      p_values[[k]][j] <- aa$p.value
    }
  }
  p_values <- do.call(cbind,p_values)
  odds_ratio <- do.call(cbind,odds_ratio)
  p_values <- apply(p_values,1,median)
  odds_ratio <- apply(odds_ratio,1,median)
  p_values <- p.adjust(p_values,method="BH")
  fisher_test_cell_types <- data.frame(cell_type=cell_types,p_values=p_values,odds_ratio=odds_ratio)
  fisher_test_cell_types <- fisher_test_cell_types[order(fisher_test_cell_types$p_values),]
  fisher_test_cell_types$sig <- "enriched"
  fisher_test_cell_types$sig[ fisher_test_cell_types$odds_ratio < 1] <- "depleted"
  fisher_test_cell_types$sig[fisher_test_cell_types$p_values > alpha] <- "not significant"

  p <- ggplot( fisher_test_cell_types , aes(x=log10(odds_ratio), y=-log10(p_values+1e-100),text=cell_type)) +
    labs(x=expression(paste("log"[10]," of odds ratio")),y=expression(paste("-log"[10]," of FDR adjusted p-values"))) +
    geom_point(aes(color=sig),size=3) +
    scale_color_manual(values=c("not significant"="grey","enriched"="darkblue","depleted" = "darkred"),name="") +
    ggrepel::geom_text_repel(data=fisher_test_cell_types[(apply(cbind(fisher_test_cell_types$odds_ratio,1/fisher_test_cell_types$odds_ratio),1,max)>1.5) | (fisher_test_cell_types$p_values < 0.1),],
                             aes(x=log10(odds_ratio), y=-log10(p_values+1e-100), label=cell_type), max.overlaps=Inf, size=4) +
    theme_classic(base_size=14) + scale_y_continuous(trans=pseudo_log_trans(sigma=1,base=10))+
    theme(
      axis.text = element_text(size=rel(0.75), color='black'),
      axis.title = element_text(size=rel(1.0), color='black')) +
    annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = -Inf, fill= "blue",alpha=0.1) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = Inf, ymax = -Inf, fill= "red",alpha=0.1)
  print(p)
  return(fisher_test_cell_types)
}



