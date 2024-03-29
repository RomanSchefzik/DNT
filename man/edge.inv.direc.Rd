% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SourceCodePTNWNeu.R
\name{edge.inv.direc}
\alias{edge.inv.direc}
\title{Calculates the differences in edge weight for each network edge}
\usage{
edge.inv.direc(A, B)
}
\arguments{
\item{A, B}{adjacency matrices of the same dimension}
}
\value{
a symmetric \eqn{N \times N} matrix, with \eqn{N} denoting the number of nodes, containing the differences in edge weight between two respective nodes
}
\description{
This function calculates the differences in edge weight for each edge between two networks (adjacency matrices of the same dimension).
}
\details{
This function calculates the differences in edge weight for each edge between two networks (adjacency matrices of the same dimension). In particular, it takes account of potentially different directions (signs) of the edge weights (associations) when deriving the differences. Differences in edge weight are one of the edge-specific network difference characteristics to compare two networks.
}
\examples{
A<-create.adjacency.matrix(ExDataA,methodlist=list("Spearman"))
B<-create.adjacency.matrix(ExDataB,methodlist=list("Spearman"))
edge.inv.direc(A,B)
}
