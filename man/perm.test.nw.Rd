% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SourceCodePTNWNeu.R
\name{perm.test.nw}
\alias{perm.test.nw}
\title{Permutation-based test for differences between two networks}
\usage{
perm.test.nw(
  A,
  B,
  permnum,
  methodlist,
  thresh = 0.05,
  score.funct,
  paired = FALSE
)
}
\arguments{
\item{A, B}{input data tables from which the adjacency matrices will be generated, to be provided in form of matrices, arrays, data frames or tibbles; need to have the same number of columns (corresponding to the number of nodes)}

\item{permnum}{a number, specifying the number of permutations}

\item{methodlist}{a list specifying the method which is used to estimate and create the adjacency matrices; see details for possible options and further information}

\item{thresh}{a number between 0 and 1 (default is set to 0.05) specifying the singificance level: if the p-value corresponding to an edge weight is greater than \code{thresh}, the corresponding edge weight is not considered to be significant and thus set to zero}

\item{score.funct}{the function used to compare the adjacency matrices A and B; see details for possible options and further information}

\item{paired}{Boolean, specifying whether the data underlying the two networks is paired or not}
}
\value{
a list, whose specific structure depends on the specified network difference characteristic in the argument \code{score.funct} (with \eqn{N} denoting the number of nodes):
\itemize{
\item{if \code{score.funct} is one of \code{frobenius.metric}, \code{global.str}, \code{maximum.metric}, \code{spec.dist} or \code{jaccard.dist}:} \cr {a list with 6 elements: the adjacency matrix for input dat set \code{A} (\eqn{N \times N} matrix), the adjacency matrix for input dat set \code{B} (\eqn{N \times N} matrix), the value of the test statistic (vector of length 1), the values of the test statistics when applying the permutations (vector of length \code{permnum}), the p-value (vector of length 1), the p-value when inserting a pseudocount in the p-value calculation to avoid p-values that are exactly zero in permutation-based settings (vector of length 1)}
\item{if \code{score.funct} is \code{number.differences}:} \cr {a list with 6 elements: output when applying \code{create.graph} to input data table \code{A} (list with 15 elements; see documentation of \code{create.graph} function for details), output when applying \code{create.graph} to input data table \code{B} (list with 15 elements; see documentation of \code{create.graph} function for details), the value of the test statistics (vector of length 6), the values of the test statistics when applying the permutations (\code{permnum}\eqn{\times 6} matrix), the p-values (vector of length 6), the p-values when inserting a pseudocount in the p-value calculation to avoid p-values that are exactly zero in permutation-based settings (vector of length 6)}
\item{if \code{score.funct} is one of \code{betweenness.inv}, \code{closeness.inv}, \code{degree.inv} or \code{eigen.inv}:} \cr {a list with 6 elements: the adjacency matrix for input dat set \code{A} (\eqn{N \times N} matrix), the adjacency matrix for input dat set \code{B} (\eqn{N \times N} matrix), the value of the test statistic for each node (vector of length \eqn{N}), the values of the test statistics for each node when applying the permutations (\code{permnum}\eqn{\times N} matrix), the p-value for each node (vector of length \eqn{N}), the p-value for each node when inserting a pseudocount in the p-value calculation to avoid p-values that are exactly zero in permutation-based settings (vector of length \eqn{N})}
\item{if \code{score.funct} is one of \code{edge.inv} or \code{edge.inv.direc}:} \cr {a list with 8 elements:  the adjacency matrix for input dat set \code{A} (\eqn{N \times N} matrix), the adjacency matrix for input dat set \code{B} (\eqn{N \times N} matrix), the value of the test statistic for each node-node pair (\eqn{N \times N} matrix), the values of the test statistics for each node-node pair when applying the permutations (\code{permnum}\eqn{\times N \times N} array), the p-value for each node-node pair (\eqn{N \times N} matrix), the p-value for each node-node pair when inserting a pseudocount in the p-value calculation to avoid p-values that are exactly zero in permutation-based settings (\eqn{N \times N} matrix), a simplified overview of the p-values for the node-node pairs (\eqn{\frac{N(N-1)}{2} \times 3} matrix; node-node pair in columns 1 and 2, corresponding p-value in column 3), a simplified overview of the p-values for the node-node pairs when inserting a pseudocount in the p-value calculation (\eqn{\frac{N(N-1)}{2} \times 3} matrix; node-node pair in columns 1 and 2, corresponding p-value in column 3)}
}
}
\description{
This function provides a permutation-based frame for testing for differences between two networks. In particular, various (i) network estimation methods and (ii) network difference characteristics can be specified.
}
\details{
This function provides a permutation-based frame for testing for differences between two networks. In particular, various (i) network estimation methods and (ii) network difference characteristics can be specified.
\cr
(i) The network estimation method has to be specified in form of a list in the \code{methodlist} argument. Currently, the following estimation methods are supported:
\itemize{
\item{\code{list("Spearman")} \cr Edge weights are estimated using Spearman correlation, where unadjusted p-values are employed to determine significance. To apply this method, the expression \code{list("Spearman")} has to be provided in the \code{methodlist} argument.}
\item{\code{list("Spearman.adj",adjustment method)} \cr Edge weights are estimated using Spearman correlation, where p-values adjusted for multiple testing are employed to determine significance. To apply this method, the expression \code{list("Spearman.adj",adjustment method)} has to be provided in the \code{methodlist} argument, where \code{adjustment method} has to be one of the options for multple testing adjustment provided by the standard \code{p.adjust} R function, i.e. one of \code{"BH"}, \code{"bonferroni"}, \code{"BY"}, \code{"fdr"}, \code{"hochberg"}, \code{"holm"} or \code{"hommel"}.}
\item{\code{list("PCSpearman")} \cr Edge weights are estimated using partial Spearman correlation, where unadjusted p-values are employed to determine significance. To apply this method, the expression \code{list("PCSpearman")} has to be provided in the \code{methodlist} argument.}
\item{\code{list("PCSpearman.adj",adjustment method)} \cr Edge weights are estimated using partial Spearman correlation, where p-values adjusted for multiple testing are employed to determine significance. To apply this method, the expression \code{list("PCSpearman.adj",adjustment method)} has to be provided in the \code{methodlist} argument, where \code{adjustment method} has to be one of the options for multple testing adjustment provided by the standard \code{p.adjust} R function, i.e. one of \code{"BH"}, \code{"bonferroni"}, \code{"BY"}, \code{"fdr"}, \code{"hochberg"}, \code{"holm"} or \code{"hommel"}.}
\item{\code{list("DistCorr")} \cr Edge weights are estimated using distance correlation, where unadjusted p-values are employed to determine significance. To apply this method, the expression \code{list("DistCorr")} has to be provided in the \code{methodlist} argument. Note that the calculations may require larger computation times, as a permutation test is involved to derive the corresponding p-values for the distance correlations.}
\item{\code{list("DistCorr.adj",adjustment method)} \cr Edge weights are estimated using distance correlation, where p-values adjusted for multiple testing are employed to determine significance. To apply this method, the expression \code{list("DistCorr.adj",adjustment method)} has to be provided in the \code{methodlist} argument, where \code{adjustment method} has to be one of the options for multple testing adjustment provided by the standard \code{p.adjust} R function, i.e. one of \code{"BH"}, \code{"bonferroni"}, \code{"BY"}, \code{"fdr"}, \code{"hochberg"}, \code{"holm"} or \code{"hommel"}. Note that the calculations may require larger computation times, as a permutation test is involved to derive the corresponding p-values for the distance correlations.}
\item{\code{list("EBICglasso",correlation type,tuning parameter)} \cr Edge weights are estimated using the EBICglasso approach. To apply this method, the expression \code{list("EBICglasso",correlation type,tuning parameter)} has to be provided in the \code{methodlist} argument. Here, \code{correlation type} has to be one of the association options provided by the standard \code{cor} R function, i.e. one of \code{"kendall"}, \code{"pearson"} or \code{"spearman"}. Moreover, \code{tuning parameter} has to be a number specifying the EBIC tuning parameter \eqn{\gamma}. Typical choices include values between 0 and 0.5, where smaller values usually lead to a higher sensitivity in that more edges are included into the network. \cr Note that for EBICglasso, an additional specification of the \code{thresh} argument is obsolete, as it is not used for the application of the method.}
}
(ii) To quantify differences between two networks, the following (a) overall, (b) edge-specific and (c) node-specific network difference characteristics, which have to be supplied in the \code{score.funct} argument, are currently supported:
\itemize{
\item{\code{frobenius.metric} (overall)} \cr {Calculates the Frobenius metric between two networks}
\item{\code{global.str} (overall)} \cr {Calculates the difference in global strength between two networks}
\item{\code{maximum.metric} (overall)} \cr {Calculates the maximum metric between two networks}
\item{\code{number.differences} (overall)} \cr {Calculates the differences in numbers of edges, clusters and isolated nodes between two networks}
\item{\code{spec.dist} (overall)} \cr {Calculates the spectral distance between two networks}
\item{\code{jaccard.dist} (overall)} \cr {Calculates the Jaccard distance between two networks}
\item{\code{betweenness.inv} (node-specific)} \cr {Calclulates the differences in betweenness between two networks for each node}
\item{\code{closeness.inv} (node-specific)} \cr {Calculates the differences in closeness between two networks for each node}
\item{\code{degree.inv} (node-specific)} \cr {Calculates the differences in degree between two networks for each node}
\item{\code{eigen.inv} (node-specific)} \cr {Calculates the differences in eigenvector centrality between two networks for each node}
\item{\code{edge.inv} (edge-specific)} \cr {Calculates the differences in absolute edge weights between two networks for each edge}
\item{\code{edge.inv.direc} (edge-specific)} \cr {Calculates the differences in edge weights between two networks for each edge}
}
Typically, a large number of permutations (e.g. 1000 or 10000) should be chosen in order to obtain reliable results. Note that a large number of permutations may lead to increasing computation times.

Note that when underlying data is paired (\code{paired=TRUE}), the input data tables need to have exactly the same dimension, i.e. the same number of columns (nodes) AND rows (samples).
}
\examples{
##examples using (a) overall network difference characteristics
res1<-perm.test.nw(A=ExDataA,B=ExDataB,permnum=10000,methodlist=list("PCSpearman"),
score.funct=frobenius.metric)
res2<-perm.test.nw(A=ExDataA,B=ExDataB,permnum=10000,methodlist=list("Spearman"),
score.funct=global.str,paired=TRUE)

##examples using (b) node-specific network difference characteristics
res3<-perm.test.nw(A=ExDataA,B=ExDataB,permnum=10000,methodlist=list("Spearman"),
thresh=0.1,score.funct=betweenness.inv,paired=TRUE)
res4<-perm.test.nw(A=ExDataA,B=ExDataB,permnum=10000,methodlist=list("EBICglasso",
"spearman",0.1), score.funct=degree.inv)

##examples using (c) edge-specific network difference characteristics
res5<-perm.test.nw(A=ExDataA,B=ExDataB,permnum=10000,methodlist=list("Spearman.adj",
"bonferroni"),score.funct=edge.inv)
res6<-perm.test.nw(A=ExDataA,B=ExDataB,permnum=10000,methodlist=list("EBICglasso",
"spearman",0.01),score.funct=edge.inv.direc,paired=TRUE)

}
