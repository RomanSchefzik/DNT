% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plots.R
\name{create.Igraphclustering}
\alias{create.Igraphclustering}
\title{Creates \code{igraph} graph and clustering results}
\usage{
create.Igraphclustering(A, methodlist, thresh = 0.05)
}
\arguments{
\item{A}{input data table from which the adjacency matrix will be generated, to be provided in form of a matrix, array, data frame or tibble}

\item{methodlist}{a list specifying the method which is used to estimate and create the adjacency matrix; see details for further information}

\item{thresh}{a number between 0 and 1 (default is set to 0.05) specifying the singificance level: if the p-value corresponding to an edge weight is greater than \code{thresh}, the corresponding edge weight is not considered to be significant and thus set to zero}
}
\value{
a list with two elements: the first list element contains the results for the clustering by the Girvan-Newman algorithm based on edge betweenness; the second list element contains the results for the creation of an \code{igraph} graph
}
\description{
This function creates \code{igraph} graph and clustering results from an adjacency matrix that is produced out of an input data table, using an estimation method specified by the user.
}
\details{
This function creates \code{igraph} graph and clustering results from an adjacency matrix that is produced out of an input data table, using an estimation method specified by the user. Clustering is performed using the Girvan-Newman algorithm based on edge betweenness. The function is used in the central \code{comp.plot} function of the DNT package for visual comparison of two networks and builds on respective implementations in the \code{igraph} R package. \cr
The network estimation method has to be specified in form of a list in the \code{methodlist} argument. Currently, the following estimation methods are supported:
\itemize{
\item{\code{list("Spearman")} \cr Edge weights are estimated using Spearman correlation, where unadjusted p-values are employed to determine significance. To apply this method, the expression \code{list("Spearman")} has to be provided in the \code{methodlist} argument.}
\item{\code{list("Spearman.adj",adjustment method)} \cr Edge weights are estimated using Spearman correlation, where p-values adjusted for multiple testing are employed to determine significance. To apply this method, the expression \code{list("Spearman.adj",adjustment method)} has to be provided in the \code{methodlist} argument, where \code{adjustment method} has to be one of the options for multple testing adjustment provided by the standard \code{p.adjust} R function, i.e. one of \code{"BH"}, \code{"bonferroni"}, \code{"BY"}, \code{"fdr"}, \code{"hochberg"}, \code{"holm"} or \code{"hommel"}.}
\item{\code{list("PCSpearman")} \cr Edge weights are estimated using partial Spearman correlation, where unadjusted p-values are employed to determine significance. To apply this method, the expression \code{list("PCSpearman")} has to be provided in the \code{methodlist} argument.}
\item{\code{list("PCSpearman.adj",adjustment method)} \cr Edge weights are estimated using partial Spearman correlation, where p-values adjusted for multiple testing are employed to determine significance. To apply this method, the expression \code{list("PCSpearman.adj",adjustment method)} has to be provided in the \code{methodlist} argument, where \code{adjustment method} has to be one of the options for multple testing adjustment provided by the standard \code{p.adjust} R function, i.e. one of \code{"BH"}, \code{"bonferroni"}, \code{"BY"}, \code{"fdr"}, \code{"hochberg"}, \code{"holm"} or \code{"hommel"}.}
\item{\code{list("DistCorr")} \cr Edge weights are estimated using distance correlation, where unadjusted p-values are employed to determine significance. To apply this method, the expression \code{list("DistCorr")} has to be provided in the \code{methodlist} argument. Note that the calculations may require larger computation times, as a permutation test is involved to derive the corresponding p-values for the distance correlations.}
\item{\code{list("DistCorr.adj",adjustment method)} \cr Edge weights are estimated using distance correlation, where p-values adjusted for multiple testing are employed to determine significance. To apply this method, the expression \code{list("DistCorr.adj",adjustment method)} has to be provided in the \code{methodlist} argument, where \code{adjustment method} has to be one of the options for multple testing adjustment provided by the standard \code{p.adjust} R function, i.e. one of \code{"BH"}, \code{"bonferroni"}, \code{"BY"}, \code{"fdr"}, \code{"hochberg"}, \code{"holm"} or \code{"hommel"}. Note that the calculations may require larger computation times, as a permutation test is involved to derive the corresponding p-values for the distance correlations.}
\item{\code{list("EBICglasso",correlation type,tuning parameter)} \cr Edge weights are estimated using the EBICglasso approach. To apply this method, the expression \code{list("EBICglasso",correlation type,tuning parameter)} has to be provided in the \code{methodlist} argument. Here, \code{correlation type} has to be one of the association options provided by the standard \code{cor} R function, i.e. one of \code{"kendall"}, \code{"pearson"} or \code{"spearman"}. Moreover, \code{tuning parameter} has to be a number specifying the EBIC tuning parameter \eqn{\gamma}. Typical choices include values between 0 and 0.5, where smaller values usually lead to a higher sensitivity in that more edges are included into the network. \cr Note that for EBICglasso, an additional specification of the \code{thresh} argument is obsolete, as it is not used for the application of the method.}
}
}
\examples{
create.Igraphclustering(ExDataA,methodlist=list("Spearman"))
create.Igraphclustering(ExDataA,methodlist=list("Spearman.adj","bonferroni"))
create.Igraphclustering(ExDataA,methodlist=list("PCSpearman"),thresh=0.1)
create.Igraphclustering(ExDataA,methodlist=list("PCSpearman.adj","BH"),thresh=0.1)
create.Igraphclustering(ExDataA,methodlist=list("DistCorr"))
create.Igraphclustering(ExDataA,methodlist=list("DistCorr.adj","bonferroni"))
create.Igraphclustering(ExDataB,methodlist=list("EBICglasso","spearman",0.1))
create.Igraphclustering(ExDataB,methodlist=list("EBICglasso","pearson",0.05))
}
