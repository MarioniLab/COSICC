% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_functions_kinetics_DE.R
\name{COSICC_kinetics}
\alias{COSICC_kinetics}
\title{COSICC_kinetics}
\usage{
COSICC_kinetics(sce_case, sce_control, alpha = 0.1)
}
\arguments{
\item{sce_case}{An object of class SingleCellExperiment for the perturbed data set
including cells for one lineage. The colData needs to include columns labelled
cell in cancer etc.) and 'pt' (pseudotime, double, a measure of progression along the
lineage trajectory.}

\item{sce_control}{An object of class SingleCellExperiment for the control data set
including cells for one lineage. The colData needs to include columns labelled
cell in cancer etc.) and 'pt' (pseudotime, double, a measure of progression along the
lineage trajectory.}

\item{alpha}{type I error cutoff for the Wilcoxon signed rank test, the default is set to 0.1}
}
\value{
a list containing
1) p_value_case: the p-value of a Wilcoxon rank-sum test determining
whether there is significant difference in developmental progression along the lineage trajectory
between marked and unmarked cells of the case data set
2)  convidence_interval_overlapping_with_control: logical, if TRUE then the 95 percent confidence intervals for
the Wilcoxon statistics for the case and control data sets overlap, and any difference between
marked and unmarked cells for the case data set cannot be seen as significant relative to the control
3)  sig: logical, whether the p-value is significant at level alpha, and there is no overlap between the
confidence intervals for the case and control data sets
}
\description{
COSICC_kinetics
}
\details{
COSICC_kinetics determines whether there is significant developmental delay
or acceleration along a lineage trajectory
}
