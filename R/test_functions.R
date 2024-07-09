
simulate_no_DE <- function(n,m){
    simulate_DE_dynamic(n,de.prob_case=0,de.prob_control=0,n_genes=m,plot=FALSE)
}

#' @title test_FDR
#' @details This function simulates data with no differential expression to test FDR control
#' for COSICC_DE.
#' @param n number of lineage trajectories to be simulated for the test data set.
#' @param m number of genes to be simulated and tested for each lineage trajectory
#' @export
test_FDR <- function(n,m){
    no_DE_simulation <- simulate_no_DE(n,m=m)
    COSICC_DE_no_DE <- COSICC_DE(no_DE_simulation$list_sce_case,no_DE_simulation$list_sce_control,
                   thresh_p_adj=0.1,
                    list_lineage_genes=rep(list(rownames(no_DE_simulation$list_sce_case[[1]])),n))
    nr_DE_static <- unlist(lapply(COSICC_DE_no_DE[1:length(no_DE_simulation$list_sce_case)],function(x) length(x$selected_static)))
    nr_DE_dynamic <- unlist(lapply(COSICC_DE_no_DE[1:length(no_DE_simulation$list_sce_case)],function(x) length(x$selected_dynamic)))
    return(list(FDR_static=nr_DE_static/(m*n),FDR_dynamic=nr_DE_dynamic/(m*n)))
}


#FDR_testing <- test_FDR(10,10)

