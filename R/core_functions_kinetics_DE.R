globalVariables(c("DC1", "DC2","celltype","dpt","effect","FDR","gene"))
#' @import scran
#' @import scater
#' @import ggthemes
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats wilcox.test
#' @importFrom stats p.adjust
#' @importFrom stats lm
#' @importFrom stats complete.cases
#' @importFrom destiny DiffusionMap
#' @importFrom destiny eigenvectors
#' @importFrom stats cor
#' @import ggplot2



#' @title find_dynamic_genes
#' @details Identification of genes that vary significantly across stages (FDR < max_FDR_across_time_points)
#' and for which there is no strong evidence of variation within batches of the same stage for a reference data set.
#' @param sce An object of class SingleCellExperiment including cells for one lineage. The colData
#' needs to include columns labelled 'stage' and 'batch'.
#' @param excluded_genes Genes to be excluded a priory from the selection of dynamic genes,
#' by default NULL (none).
#' @param FDR_across_time_points see above
#' @return An object of type SingleCellExperiment subset to the dynamic genes identified, and otherwise
#' identical to the input sce.
#' @export
find_dynamic_genes <- function(sce,excluded_genes=NULL,FDR_across_time_points=10^-10)
{
    sce <- sce[,!(is.na(sce$stage))]
    stage_vec <- sce$stage
    lf_init_end <- rowMeans(logcounts(sce)[,stage_vec == max(stage_vec)])- rowMeans(logcounts(sce)[,stage_vec == min(stage_vec)])
    sce <- sce[abs(lf_init_end) > 0.5,]
    stage_vec <- sce$stage
    p_vals <- p.adjust(apply(logcounts(sce),1,function(v) summary(lm(v ~ stage_vec))$coefficients[2,4]),method="BH")
    residuals <- apply(logcounts(sce),1,function(v) lm(v ~ stage_vec)$residuals)
    sce <- sce[p_vals<FDR_across_time_points,]
    residuals <- residuals[,p_vals<FDR_across_time_points]
    p_vals_res <- p.adjust(apply(residuals,2,function(v) summary(lm(v ~ sce$batch))$coefficients[2,4]),method="BH")
    sce <- sce[p_vals_res > 0.2,]
    return(sce)
}



#' @title compute_and_plot_pseudotime
#' @details The function compute_and_plot_pseudotime computes diffusion maps, and identifies the diffusion component
#' most correlated with embryonic stage. The function also prints plots of the first two diffusion components,
#' @param sce An object of type SingleCellExperiment including cells for one lineage, the colData of each of the SingleCellExperiments
#' in the list needs to include the slot 'stage' for embryonic stage.
#' @export
compute_and_plot_pseudotime <- function(sce){
    dm <- DiffusionMap(t(as.matrix(logcounts(sce))),n_pcs=10)
    cors <- cor(dm@eigenvectors,sce$stage)
    dm@eigenvectors[,cors <0 ] <- -dm@eigenvectors[,cors <0 ]
    sce$dpt <- dm@eigenvectors[,which.max(abs(cors))]
    xx <- sample(1:length(dm$DC1))
    tmp <- data.frame(DC1 = eigenvectors(dm)[xx, 1],
                      DC2 = eigenvectors(dm)[xx, 2],
                      celltype = sce$celltype[xx],
                      dpt = sce$dpt[xx])
    p1 <- ggplot(tmp, aes(x = DC1, y = DC2, colour = celltype)) +
        geom_point() +
        xlab("Diffusion component 1") +
        ylab("Diffusion component 2") +
        theme_classic(base_size=11) + theme(legend.position = "bottom",legend.box="vertical",legend.margin=margin())+
        labs(color="")+ guides(color = guide_legend(nrow=2,override.aes = list(size = 3)))+ scale_color_colorblind()
    print(p1)

    p2 <- ggplot(tmp, aes(x = DC1, y = DC2, colour = dpt)) +
        geom_point(alpha=0.5,size=2) +
        xlab("Diffusion component 1") +
        ylab("Diffusion component 2") +
        theme_classic(base_size=11) +scale_color_viridis_c()+ theme(legend.position = "bottom",legend.box="vertical",legend.margin=margin())+
        labs(color="")+ guides(color = guide_legend(override.aes = list(size = 3)))
    print(p2)
    return(sce)
}


correlation_pseudotime <- function(reference_sce,perturbed_sce){
    correlation_matrix <- matrix(ncol=ncol(reference_sce),nrow=ncol(perturbed_sce))
    rownames(correlation_matrix) <- colnames(perturbed_sce)
    colnames(correlation_matrix) <- colnames(reference_sce)
    logcounts_reference <- logcounts(reference_sce)
    logcounts_perturbed <- logcounts(perturbed_sce)
    for (j in 1:ncol(correlation_matrix)){
        for (k in 1:nrow(correlation_matrix)){
            correlation_matrix[k,j] <- cor(logcounts_reference[,j],logcounts_perturbed[,k])
        }
    }
    perturbed_sce$dpt <- apply(correlation_matrix, 1, function(x)
        mean(reference_sce$dpt[colnames(correlation_matrix)[order(x,decreasing=T) %in% 1:10]]))
    output <- c()
    output$correlation_matrix <- correlation_matrix
    return(output)
}

Wilcoxon_test_marked_vs_unmarked <- function(marked_sce){
    wilcox_test <- wilcox.test(marked_sce$pt[marked_sce$marked],
                               marked_sce$pt[!(marked_sce$marked)],conf.int = TRUE)
    return(data.frame(lower_bound = wilcox_test$conf.int[1],
                      upper_bound = wilcox_test$conf.int[2],
                      p_value = wilcox_test$p.value,
                      estimate = wilcox_test$estimate ))
}


#' @title COSICC_kinetics
#' @details COSICC_kinetics determines whether there is significant developmental delay
#' or acceleration along a lineage trajectory
#' @param sce_case An object of class SingleCellExperiment for the perturbed data set
#' including cells for one lineage. The colData needs to include columns labelled
# 'marked' (logical, indicating whether the cell is marked by a fluorescent marker, or a malignant
#' cell in cancer etc.) and 'pt' (pseudotime, double, a measure of progression along the
#' lineage trajectory.
#' @param sce_control An object of class SingleCellExperiment for the control data set
#' including cells for one lineage. The colData needs to include columns labelled
# 'marked' (logical, indicating whether the cell is marked by a fluorescent marker, or a malignant
#' cell in cancer etc.) and 'pt' (pseudotime, double, a measure of progression along the
#' lineage trajectory.
#' @param alpha type I error cutoff for the Wilcoxon signed rank test, the default is set to 0.1
#' @return a list containing
#' 1) p_value: the p-value of a Wilcoxon rank-sum test determining
#' whether there is significant difference in developmental progression along the lineage trajectory
#' between marked and unmarked cells of the case data set
#' 2)  convidence_interval_overlapping_with_control: logical, if TRUE then the 95 percent confidence intervals for
#' the Wilcoxon statistics for the case and control data sets overlap, and any difference between
#' marked and unmarked cells for the case data set cannot be seen as significant relative to the control
#' 3)  sig: logical, whether the p-value is significant at level alpha, and there is no overlap between the
#' confidence intervals for the case and control data sets
#' @export
COSICC_kinetics <- function(sce_case,sce_control,alpha=0.1){
    Wilcoxon_test_case <- Wilcoxon_test_marked_vs_unmarked(sce_case)
    Wilcoxon_test_control <- Wilcoxon_test_marked_vs_unmarked(sce_control)
    overlapping_with_control <- !(Wilcoxon_test_case$lower_bound > Wilcoxon_test_control$upper_bound |
                                      Wilcoxon_test_case$upper_bound < Wilcoxon_test_control$lower_bound)
    out <- list(p_value=Wilcoxon_test_case$p_value,
                convidence_interval_overlapping_with_control =
                    overlapping_with_control,
                sig = Wilcoxon_test_case$p_value < alpha & !(overlapping_with_control))
    return(out)
}


volcano_plot <- function(vector_FDR,vector_effect_sizes,thresh_FDR=0.1,max_highlight=20,highlight_extra_genes=NULL,
                         x_label="coefficient",title=""){
    vector_effect_sizes <- vector_effect_sizes[names(vector_FDR)]
    df <- data.frame(FDR=vector_FDR,effect=vector_effect_sizes)
    df <- df[complete.cases(df),]
    df$gene <- rownames(df)
    df$sig <- "not significant"
    df$sig[df$FDR<thresh_FDR] <- "significant"
    max_highlight <- min(nrow(df),max_highlight)
    highlight <- df$FDR < thresh_FDR &
        ((df$effect > 0 & df$effect >= sort(df$effect,decreasing=TRUE)[max_highlight])|(df$effect < 0 & df$effect <= sort(df$effect,decreasing=FALSE)[max_highlight]))
    highlight[highlight & (df$effect > 0)] <- "up"
    highlight[highlight==TRUE & (df$effect < 0)] <- "down"
    highlight[df$gene %in% highlight_extra_genes] <- "gene_of_interest"
    df$highlight <- highlight
    pp <- ggplot(df, aes(x=effect, y=-log10(FDR),color=highlight)) +
        geom_point(size=3, alpha=0.4) +
        ggtitle(title) + scale_color_manual(values=c("FALSE"="grey","up"="darkblue","down"="darkred","gene_of_interest"="purple"))+
        labs(y=expression('-Log'[10]*' FDR'), x=x_label) +
        theme_classic(base_size=20) +
        theme(legend.position="none", plot.title = element_text(size = rel(1), hjust = 0.5))+
        geom_text_repel(data=df[!(highlight==FALSE),],
                        aes(x = effect, y = -log10(FDR),label=gene),max.overlaps=100)+
        scale_y_continuous(trans=scales::pseudo_log_trans(sigma=5,base = 1.5))+
        geom_hline(yintercept = 1)
    return(pp)

}


