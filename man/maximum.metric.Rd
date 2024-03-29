% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SourceCodePTNWNeu.R
\name{maximum.metric}
\alias{maximum.metric}
\title{Calculates the maximum metric}
\usage{
maximum.metric(A, B)
}
\arguments{
\item{A, B}{adjacency matrices of the same dimension}
}
\value{
the maximum metric between the two adjacency matrices \code{A} and \code{B}
}
\description{
This function calculates the maximum metric between two adjacency matrices of the same dimension.
}
\details{
This function calculates the maximum metric between two adjacency matrices of the same dimension. The maximum metric is an overall characteristic that can be employed to compare two networks.
}
\examples{
A<-create.adjacency.matrix(ExDataA,methodlist=list("Spearman"))
B<-create.adjacency.matrix(ExDataB,methodlist=list("Spearman"))
maximum.metric(A,B)
}
