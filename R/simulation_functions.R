#' @import splatter
simulate_DE_dynamic_core <- function(seed,de.prob_case,de.prob_control,n_genes){
    sce_DE_dynamic_control_marked <- splatSimulatePaths(batchCells=250,verbose=FALSE,nGenes=n_genes,
                        de.prob=de.prob_control, seed=seed,de.facLoc=0.2,
                           path.nSteps=10,path.nonlinearProb=0)
    sce_DE_dynamic_control_unmarked <-
        splatSimulatePaths(batchCells=250,verbose=FALSE,nGenes=n_genes,de.prob=0,
                           seed=seed,path.nSteps=10,path.nonlinearProb=0)
    rowData(sce_DE_dynamic_control_marked)[,5:7] <- NULL
    rowData(sce_DE_dynamic_control_unmarked)[,5:7] <- NULL
    sce_DE_dynamic_control <- cbind(sce_DE_dynamic_control_unmarked,sce_DE_dynamic_control_marked)
    colnames(sce_DE_dynamic_control) <- paste0("Gene",1:ncol(sce_DE_dynamic_control))

    sce_DE_dynamic_control$pt <- as.double(sce_DE_dynamic_control$Step)
    sce_DE_dynamic_control$Step <- NULL
    sce_DE_dynamic_control <- logNormCounts(sce_DE_dynamic_control)

    sce_DE_dynamic_control$marked <- c(rep(FALSE,ncol(sce_DE_dynamic_control_unmarked)),
                                       rep(TRUE,ncol(sce_DE_dynamic_control_marked)))
    sce_DE_dynamic_case_marked <-
        splatSimulatePaths(batchCells=250,verbose=FALSE,nGenes=n_genes,de.prob=de.prob_case, seed=seed,de.facLoc=0.6,
                           path.nSteps=10,path.nonlinearProb=0)
    sce_DE_dynamic_case_unmarked <-
        splatSimulatePaths(batchCells=250,verbose=FALSE,nGenes=n_genes,de.prob=0, seed=seed,
                         de.facLoc=0.4,path.nSteps=10,path.nonlinearProb=0)
    rowData(sce_DE_dynamic_case_marked)[,5:7] <- NULL
    rowData(sce_DE_dynamic_case_unmarked)[,5:7] <- NULL
    sce_DE_dynamic_case <- cbind(sce_DE_dynamic_case_unmarked,sce_DE_dynamic_case_marked)
    colnames(sce_DE_dynamic_case) <- paste0("Gene",1:ncol(sce_DE_dynamic_case))

    sce_DE_dynamic_case$pt <- as.double(sce_DE_dynamic_case$Step)
    sce_DE_dynamic_case$Step <- NULL
    sce_DE_dynamic_case <- logNormCounts(sce_DE_dynamic_case)

    sce_DE_dynamic_case$marked <- c(rep(FALSE,ncol(sce_DE_dynamic_case_unmarked)),
                                    rep(TRUE,ncol(sce_DE_dynamic_case_marked)))

    return(list(sce_case=sce_DE_dynamic_case,sce_control=sce_DE_dynamic_control))
}

#' @import gridExtra
#' @title simulate_DE_dynamic
#' @description Simulation of dynamic differential expression, i.e differences between marked and
#' unmarked cells in terms of changes in expression along a lineage trajectory.
#' @param n integer, number of lineages to be simulated
#' @return A nested list consisting of the following three lists
#' @param list_sce_case list of SingleCellExperiments for different lineages
#' of a simulated perturbation experiment
#' @param list_sce_control list of corresponding control experiments
#' @param de.prob_case probability of differential expression for marked cells
#' compared to unmarked cells for the case data set
#' @param de.prob_control probability of differential expression for marked cells
#' compared to unmarked cells for the control data set
#' @param n_genes number of genes to be simulated
#' @param plot whether to produce gene expression plots for simulated genes
#' @return A list containing the following two lists
#' @param list_sce_case list of simulated perturbation data sets
#' @param list_sce_control list of simulated control data sets
#' @export
simulate_DE_dynamic <- function(n,seeds=NULL,de.prob_case=0.35,de.prob_control=0.2,n_genes=4,
                                plot=TRUE){
    if (is.null(seeds)){
        seeds <- sample(1:10000,n)
    }
    temp <- lapply(1:n,function(x) simulate_DE_dynamic_core(seeds[x],de.prob_case=de.prob_case,
                de.prob_control=de.prob_control,n_genes=n_genes))
    list_sce_case <- lapply(temp,function(x) return (x$sce_case))
    names(list_sce_case) <- paste0("lineage_trajectory_",1:n)
    list_sce_control <- lapply(temp,function(x) return (x$sce_control))
    names(list_sce_control) <- paste0("lineage_trajectory_",1:n)
    for (j in 1:length(list_sce_control)){
        plot_list_control <- list()
        plot_list_case <- list()

        for (k in 1:nrow(list_sce_case[[1]])){
            df_DE_dynamic_control <- data.frame(logcounts=logcounts(list_sce_control[[j]])[k,],
                                                marked = list_sce_control[[j]]$marked,
                                                pseudotime = as.factor(list_sce_control[[j]]$pt))
            if (plot==TRUE){
                plot_list_control[[k]] <- ggplot(df_DE_dynamic_control,
                                                 mapping=aes(x=pseudotime,y=logcounts,color=marked)) +
                    geom_jitter() +
                    scale_color_manual(values=c("FALSE"="darkblue","TRUE"="red"))+
                    ggtitle(paste0("Control: lineage trajectory ",j,", gene",k)) + theme_classic()
                df_DE_dynamic_case <- data.frame(logcounts=logcounts(list_sce_case[[j]])[k,],
                                                 marked = list_sce_case[[j]]$marked,
                                                 pseudotime = as.factor(list_sce_case[[j]]$pt))
                plot_list_case[[k]] <- ggplot(df_DE_dynamic_case,mapping=aes(x=pseudotime,y=logcounts,color=marked)) +
                    geom_jitter()+
                    scale_color_manual(values=c("FALSE"="darkblue","TRUE"="red"))+
                    ggtitle(paste0("Case: lineage trajectory ",j,", gene",k)) + theme_classic()


            }
        }
        if (plot == TRUE){
            nCol <- 2
            do.call("grid.arrange", c(plot_list_case, ncol=nCol))
            do.call("grid.arrange", c(plot_list_control, ncol=nCol))
        }

    }
    return(list(list_sce_case=list_sce_case,list_sce_control=list_sce_control))
}
