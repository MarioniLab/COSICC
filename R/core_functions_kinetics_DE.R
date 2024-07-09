#' @import SingleCellExperiment
#' @import scran
#' @import scater
#' @import ggrepel
#' @import ggthemes

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


#' @import destiny
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
#' @param p_value the p-value of a Wilcoxon rank-sum test determining
#' whether there is significant difference in developmental progression along the lineage trajectory
#' between marked and unmarked cells of the case data set
#' @param convidence_interval_overlapping_with_control logical, if TRUE then the 95 percent confidence intervals for
#' the Wilcoxon statistics for the case and control data sets overlap, and any difference between
#' marked and unmarked cells for the case data set cannot be seen as significant relative to the control
#' @param sig logical, whether the p-value is significant at level alpha, and there is no overlap between the
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


.adjusted_nonaggregated <- function(df,size_factors,dispersion){
    # following Meinshausen et al. (2009) JASA
    model_matrix <- model.matrix(gex ~ dataset_case + marked + pseudotime +
                                     dataset_case:marked + dataset_case:pseudotime +
                                     marked:pseudotime + dataset_case:marked:pseudotime, data=df)
    coefficient_name_dynamic <-"dataset_caseTRUE:markedTRUE:pseudotime"
    coefficient_name_static <- "dataset_caseTRUE:markedTRUE"
    p_val_selection <- rep(1,2)
    names(p_val_selection) <- c(coefficient_name_static,coefficient_name_dynamic)
    effect_size_temp <- rep(NA,2)
    names(effect_size_temp) <- c(coefficient_name_static,coefficient_name_dynamic)
    df_model <- as.data.frame(model_matrix[,-1])
    df_model$gex <- df$gex
    xx <- c(sample(which(df_model$dataset_caseTRUE & df_model$markedTRUE),
                   floor(0.5*sum(df_model$dataset_caseTRUE & df_model$markedTRUE)+runif(1))),
            sample(which(!(df_model$dataset_caseTRUE) & df_model$markedTRUE),
                   floor(0.5*sum(!(df_model$dataset_caseTRUE) & df_model$markedTRUE)+runif(1))),
            sample(which(df_model$dataset_caseTRUE & !(df_model$markedTRUE)),
                   floor(0.5*sum(df_model$dataset_caseTRUE & !(df_model$markedTRUE))+runif(1))),
            sample(which(!(df_model$dataset_caseTRUE) & !(df_model$markedTRUE)),
                   floor(0.5*sum(!(df_model$dataset_caseTRUE) & !(df_model$markedTRUE))+runif(1))))
    #xx <- sample(c(TRUE,FALSE),nrow(df),replace=TRUE)
    df_select_coeff <- df_model[xx,]#50% of the cells for
    #selection of predictors
    df_p_val <- df_model[-xx,]#50% of the cells to compute the p-values
    if(sd(df_select_coeff$gex) > 0.1){
        fit <- glmreg(gex ~ ., data=df_select_coeff,alpha=1,
                      nlambda=20,family="negbin",theta=1/dispersion,offset=size_factors[xx])

        coeff_lowest_bic <- fit$beta[,which.min(fit$bic),drop=FALSE]
        coeff_included <- rownames(coeff_lowest_bic)[coeff_lowest_bic[,1]!=0]
        coeff_included <- sapply(coeff_included,function(x) gsub("`","",x))

        effect_size_temp <- rep(NA,2)
        names(effect_size_temp) <- c(coefficient_name_static,coefficient_name_dynamic)
        p_val_temp <- rep(NA,2)
        names(p_val_temp) <- c(coefficient_name_static,coefficient_name_dynamic)

        if (length(coeff_included) > 0){
            df <- df_p_val[,c("gex",coeff_included),drop=FALSE]
            fit2 <- glm(gex ~  .,data=df,family=negative.binomial(1/dispersion),offset=size_factors[-xx])
            coeff_2 <- summary(fit2)$coefficients
            rownames(coeff_2) <- sapply(rownames(coeff_2),function(x) gsub("`","",x))
            coeff_temp <- coeff_2[intersect(names(effect_size_temp),rownames(coeff_2)),,drop=FALSE]
            effect_size_temp[intersect(names(effect_size_temp),rownames(coeff_2))] <- coeff_temp[,1]
            p_val_temp[intersect(names(effect_size_temp),rownames(coeff_2))] <- coeff_temp[,4]
        }
    nr_selected <- sum(!(is.na(p_val_temp)))
    p_val_temp[is.na(p_val_temp)] <- 1
    if (nr_selected != 0){
        p_val_selection <- sapply(p_val_temp,function(x) min(1,nr_selected*x))
    }else{
        p_val_selection <- p_val_temp
    }
    }
    return(list(p_values=p_val_selection,effect_size=effect_size_temp))
}

.aggregate_p_vals <- function(p_vals){
    gamma = seq(0.05, 1 - 1 / length(p_vals), by = 1/length(p_vals))
    l <- length(p_vals)
    penalty <- 1 - log(min(gamma))
    p_vals_agg <- rep(1,2)
    for(j in 1:2) {
        gamma_quantiles <- sapply(gamma,function(x) quantile(p_vals[,j]/x, x))
        gamma_quantiles <- sapply(gamma_quantiles, function(x) min(1,x))
        p_vals_agg[j] <- min(gamma_quantiles * penalty, 1)
    }
    return(p_vals_agg)
}

.rowMedians <- function(M){
    output <- c(NA,NA)
    for (j in 1:ncol(M)){
        if (any(!(is.na(M[,j])))){
            output[j] <- median(M[,j],na.rm=TRUE)
        }
    }
    return(output)
}

#' @import DESeq2
.estimate_dispersions <- function(sce){
    counts(sce) <- as.matrix(counts(sce))
    assays(sce)[[1]] <- counts(sce)
    names(assays(sce))[1] <- "counts"
    dds <- DESeqDataSet(sce, design = ~ 1)
    return(rowData(DESeq2::estimateDispersions(dds,quiet=TRUE))$dispersion)
}



#' @import MASS
#' @import mpath
#' @import batchelor
#' @title COSICC_DE_core
#' @details COSICC_DE_core uses lasso regression to determine whether a gene
#' is statically or dynamically dynamically expressed (DE) along a lineage trajectory.
#' p-values are computed following Meinshausen et al. (2009) JASA. By static differential
#' expresssion we mean that the gene is DE between marked and
#' unmarked cells, independently of the temporal or pseudotemporal state of the
#' specific cell. Dynamic DE, by contrast, refers to genes whose expression changes as a function of
#' (pseudo)time in different ways between marked and unmarked cells.
#' @param sce_case Object of class SingleCellExperiment for the perturbation data set
#' including cells for one lineage trajectory. The colData needs to include the
#' following columns: pt: observed time or pseudotime; marked: whether the
#' cell has a fluorescent or other marker (TRUE) or not (FALSE)
#' @param sce_case Object of class SingleCellExperiment for the control data set
#' including cells for one lineage trajectory. The colData needs to include the
#' following column: pt: observed time or pseudotime; marked: whether the
#' cell has a fluorescent or other marker (TRUE) or not (FALSE)
#' @param lineage_genes genes to be tested, externally determined genes that are important
#' to the progression along the lineage trajectory
#' @param N number of times the data set is split into two halves, where one half
#' is used for the identification of predictors using lasso, and the other half is
#' used for the computation of p-values. This avoid randomness in the p-values.
#' @return A list containing the following:
#' @param effect_size matrix with columns static and dynamic, one row for each gene.
#' The first column lists the log-fold-changes for static differential expression,
#' the second the effect sizes for the dynamic case.
#' @param p_val_adj_static adjusted p-values for static DE, based on correction
#' across all genes for the given lineage
#' @param p_val_adj_dynamic adjusted p-values for dynamic DE, based on correction
#' across all genes for the given lineage
#' @export
COSICC_DE_core <- function(sce_case,sce_control,lineage_genes,N=30){
    pseudotime_case <- sce_case$pt
    names(pseudotime_case) <- colnames(sce_case)
    pseudotime_control <- sce_control$pt
    names(pseudotime_control) <- colnames(sce_control)
    dataset <- c(rep("case",ncol(sce_case)),rep("control",ncol(sce_control)))
    dataset_case <- dataset == "case"
    marked <- c(sce_case$marked,sce_control$marked)
    pseudotime <- c(pseudotime_case,pseudotime_control)
    nr_time_groups <- length(unique(pseudotime))
    sce_case <- computeSumFactors(sce_case)
    sce_control <- computeSumFactors(sce_control)
    sce_list <- multiBatchNorm(list(sce_case,sce_control))
    sce_control <- sce_list[[2]]
    sce_case <- sce_list[[1]]
    sizeFactors <- c(sce_case$sizeFactor,sce_control$sizeFactor)
    disp <- .estimate_dispersions(sce_control)
    names(disp) <- rownames(sce_control)

    disp <- disp[lineage_genes]
    sce_case <- sce_case[lineage_genes,]
    sce_control <- sce_control[lineage_genes,]


    p_val_selection_static <- list()
    p_val_selection_dynamic <- list()
    effect_size  <- matrix(data=rep(NA,nrow(sce_case)*2),nrow=nrow(sce_case),
                           ncol=2)
    rownames(effect_size) <- rownames(sce_case)
    colnames(effect_size) <- c("static","dynamic")
    counts <- as.matrix(cbind(counts(sce_case),counts(sce_control)))


    counts <- counts[apply(counts,1,function(x) sum(x>0) >= 10),]

    for (j in 1:nrow(counts)){

        gex <- counts[j,]
        df <- data.frame(gex=gex,dataset_case=dataset_case,marked=marked,
                         pseudotime=pseudotime)

        temp <- lapply(1:N,function(N) .adjusted_nonaggregated(df,sizeFactors,disp[j]))
        temp_p_vals <- cbind(sapply(temp,function(x) x$p_values[1]),
                             sapply(temp,function(x) x$p_values[2]))
        rownames(temp_p_vals) <- NULL
        colnames(temp_p_vals) <- c("static","dynamic")
        aggregated_p_vals <- .aggregate_p_vals(temp_p_vals)
        p_val_selection_static[[j]] <- aggregated_p_vals[1]
        p_val_selection_dynamic[[j]] <- aggregated_p_vals[2]
        temp_lfc <- cbind(sapply(temp,function(x) x$effect_size[1]),
                          sapply(temp,function(x) x$effect_size[2]))
        rownames(temp_lfc) <- NULL
        colnames(temp_lfc) <- c("static","dynamic")

        effect_size[j,] <- .rowMedians(temp_lfc)
    }
    names(p_val_selection_static) <- rownames(sce_case)
    names(p_val_selection_dynamic) <- rownames(sce_case)
    nr_chosen <- sum(p_val_selection_dynamic < 0.1) + sum(p_val_selection_static < 0.1)
    if (nr_chosen == 0){
        p_val_adj_static <- rep(NA,length(p_val_selection_static))
        p_val_adj_dynamic <- rep(NA,length(p_val_selection_dynamic))
    }else{
        p_val_all_adj <- p.adjust(c(unlist(p_val_selection_static),unlist(p_val_selection_dynamic)),method="BY")
        p_val_adj_static <- p_val_all_adj[1:length(p_val_selection_static)]
        p_val_adj_dynamic <- p_val_all_adj[-(1:length(p_val_selection_static))]
        # p_val_adj_static <- unlist(p_val_selection_static)*nr_chosen
        # p_val_adj_static <- sapply(p_val_adj_static,function(x) min(x,1))
        # p_val_adj_dynamic <- unlist(p_val_selection_dynamic) * nr_chosen
        # p_val_adj_dynamic <- sapply(p_val_adj_dynamic,function(x) min(x,1))
    }
    output <- list(effect_size=effect_size,p_value_static=unlist(p_val_selection_static),
                   p_value_dynamic=unlist(p_val_selection_dynamic),
                   p_val_adj_static=p_val_adj_static, p_val_adj_dynamic= p_val_adj_dynamic)
    return(output)
}


#' @title combine_lineages_COSICC_DE
#' @details This function combines and plots outputs from COSICC_DE_core across differential lineage trajectories.
#' @param list_of_COSICC_DE_core_outputs list of outputs from COSICC_DE_core, one list element per lineage
#' @param names_lineage_trajectories names of the lineage trajectories associated with the output in
#' list_of_COSICC_DE_core_outputs
#' @param thresh_p_adj threshold for the adjusted p-value for selecting a gene as DE for a lineage trajectory
#' @return input list list_of_COSICC_DE_core_outputs, barplots of numbers of DE genes for each lineage for static DE,
#' (p_static) and dynamic DE (p_dynamic), list of volcano plots for static DE (volcano_plot_list_static)
#' and dynamic DE (volcano_plot_list_dynamic)
#' @export
combine_lineages_COSICC_DE <- function(list_of_COSICC_DE_core_outputs,names_lineage_trajectories,thresh_p_adj=0.1){
    names(list_of_COSICC_DE_core_outputs) <- names_lineage_trajectories
    list_of_COSICC_DE_core_outputs <- lapply(list_of_COSICC_DE_core_outputs,function(x)
    {x$selected <- x$p_value_adj_static < thresh_p_adj | x$p_value_adj_dynamic < thresh_p_adj; return(x)})
    list_of_COSICC_DE_core_outputs <- lapply(list_of_COSICC_DE_core_outputs,function(x)
    {x$selected_static <- x$p_value_adj_static < thresh_p_adj ; return(x)})
    list_of_COSICC_DE_core_outputs <- lapply(list_of_COSICC_DE_core_outputs,function(x)
    {x$selected_dynamic <- x$p_value_adj_dynamic < thresh_p_adj ; return(x)})
    # plot nr of DE genes per lineage trajectory
    df_nr_DE_genes <- data.frame(nr_DE_genes_static <- sapply(list_of_COSICC_DE_core_outputs,function(x) sum(x$p_val_adj_static < thresh_p_adj)),
                                 nr_DE_genes_dynamic <- sapply(list_of_COSICC_DE_core_outputs,function(x) sum(x$p_val_adj_dynamic < thresh_p_adj)),
                                 lineage_trajectory <- names(list_of_COSICC_DE_core_outputs))
    colnames(df_nr_DE_genes) <- c("nr_DE_genes_static","nr_DE_genes_dynamic","lineage_trajectory")
    list_of_COSICC_DE_core_outputs$p_static <- ggplot(df_nr_DE_genes[df_nr_DE_genes$nr_DE_genes_static>0,],aes(x=lineage_trajectory,y=nr_DE_genes_static)) +
        geom_bar(stat="identity") + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        xlab("") + ylab("number of genes") + ggtitle("static")
    print(list_of_COSICC_DE_core_outputs$p_static)
    list_of_COSICC_DE_core_outputs$p_dynamic <- ggplot(df_nr_DE_genes[df_nr_DE_genes$nr_DE_genes_dynamic>0,],aes(x=lineage_trajectory,y=nr_DE_genes_dynamic)) +
        geom_bar(stat="identity") + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        xlab("") + ylab("number of genes") + ggtitle("dynamic")
    print(list_of_COSICC_DE_core_outputs$p_dynamic)
    # volcano plots for DE
    list_of_COSICC_DE_core_outputs$volcano_plot_list_static <- list()
    list_of_COSICC_DE_core_outputs$volcano_plot_list_dynamic <- list()
    any_selected_dynamic <-
        names(list_of_COSICC_DE_core_outputs)[sapply(list_of_COSICC_DE_core_outputs,function(x) any(x$p_value_dynamic < thresh_p_adj))]
    any_selected_static <-
        names(list_of_COSICC_DE_core_outputs)[sapply(list_of_COSICC_DE_core_outputs,function(x) any(x$p_value_static < thresh_p_adj))]
    if (length(any_selected_static) >= 1){
        for (j in 1:length(any_selected_static)){
            list_of_COSICC_DE_core_outputs$volcano_plot_list_static[[j]] <-
                volcano_plot(vector_FDR=list_of_COSICC_DE_core_outputs[[any_selected_static[j]]]$p_val_adj_static,
                             vector_effect_sizes=list_of_COSICC_DE_core_outputs[[any_selected_static[j]]]$effect_size[,"static"],
                             thresh_FDR=0.1,max_highlight=40,highlight_extra_genes=NULL,
                             x_label="log-fold change",title=paste0(any_selected_static[j],": static"))
            print(list_of_COSICC_DE_core_outputs$volcano_plot_list_static[[j]] )
        }
    }
    if (length(any_selected_dynamic) >= 1){
        for (j in 1:length(any_selected_dynamic)){
            list_of_COSICC_DE_core_outputs$volcano_plot_list_dynamic[[j]] <-
                volcano_plot(vector_FDR=list_of_COSICC_DE_core_outputs[[any_selected_dynamic[j]]]$p_val_adj_dynamic,
                             vector_effect_sizes=list_of_COSICC_DE_core_outputs[[any_selected_dynamic[j]]]$effect_size[,"dynamic"],
                             thresh_FDR=0.1,max_highlight=20,highlight_extra_genes=NULL,
                             x_label="coefficient",title=paste0(any_selected_dynamic[j],": dynamic"))
            print(list_of_COSICC_DE_core_outputs$volcano_plot_list_dynamic[[j]])
        }
    }

    return(list_of_COSICC_DE_core_outputs)
}

#' @title COSICC_DE
#' @title COSICC_DE integrates the functions COSICC_DE_core and combine_lineages_COSICC_DE.
#' @param list_sce_case A list of objects of class SingleCellExperiment for the perturbation data set
#' including cells for one lineage trajectory, containing one such SingleCellExperiment object for each lineage.
#' The names of the list are the names of the lineages.
#' The colData for each of the SingleCellExperiment objects needs to include columns labelled
# 'marked' (logical, indicating whether the cell is marked by a fluorescent marker, or a malignant
#' cell in cancer etc.), 'lineage' and 'pt' (pseudotime, double, a measure of progression along the
#' lineage trajectory.
#' @param list_sce_control Like list_sce_case, for the control instead of the perturbation data set
#' @param thresh_p_adj double between 0 and 0.1, significance threshold for the adjusted p-values,
#' for selecting a gene as DE for a lineage trajectory
#' @param list_lineage_genes list of genes to be tested (ordered in the same way as list_sce_case and list_sce_control),
#' externally determined genes that are important to the progression along the lineage trajectories
#' @return see combine_lineages_COSICC_DE
#' @export
COSICC_DE <- function(list_sce_case,list_sce_control,thresh_p_adj=0.1,list_lineage_genes){
    list_of_outputs <- c()
    for (j in 1:length(list_sce_case)){
        list_of_outputs[[j]] <- COSICC_DE_core(sce_case=list_sce_case[[j]],
                                               sce_control=list_sce_control[[j]],list_lineage_genes[[j]])
    }
    combine_lineages_COSICC_DE(list_of_outputs,names(list_sce_case),thresh_p_adj=thresh_p_adj)
}


