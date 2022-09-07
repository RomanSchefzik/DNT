###1.) Network estimation

#' Creates adjacency matrix
#' @description
#' This function creates an adjacency matrix out of an input data table, using an estimation method specified by the user.
#' @param A input data table from which the adjacency matrix will be generated, to be provided in form of a matrix, array, data frame or tibble
#' @param methodlist a list specifying the method which is used to estimate and create the adjacency matrix; see details for further information
#' @param thresh a number between 0 and 1 (default is set to 0.05) specifying the singificance level: if the p-value corresponding to an edge weight is greater than \code{thresh}, the corresponding edge weight is not considered to be significant and thus set to zero
#' @details
#' This function creates an adjacency matrix out of an input data table, using an estimation method specified by the user. The network estimation method has to be specified in form of a list in the \code{methodlist} argument. Currently, the following estimation methods are supported:
#' \itemize{
#' \item{\code{list("Spearman")} \cr Edge weights are estimated using Spearman correlation, where unadjusted p-values are employed to determine significance. To apply this method, the expression \code{list("Spearman")} has to be provided in the \code{methodlist} argument.}
#' \item{\code{list("Spearman.adj",adjustment method)} \cr Edge weights are estimated using Spearman correlation, where p-values adjusted for multiple testing are employed to determine significance. To apply this method, the expression \code{list("Spearman.adj",adjustment method)} has to be provided in the \code{methodlist} argument, where \code{adjustment method} has to be one of the options for multple testing adjustment provided by the standard \code{p.adjust} R function, i.e. one of \code{"BH"}, \code{"bonferroni"}, \code{"BY"}, \code{"fdr"}, \code{"hochberg"}, \code{"holm"} or \code{"hommel"}.}
#' \item{\code{list("PCSpearman")} \cr Edge weights are estimated using partial Spearman correlation, where unadjusted p-values are employed to determine significance. To apply this method, the expression \code{list("PCSpearman")} has to be provided in the \code{methodlist} argument.}
#' \item{\code{list("PCSpearman.adj",adjustment method)} \cr Edge weights are estimated using partial Spearman correlation, where p-values adjusted for multiple testing are employed to determine significance. To apply this method, the expression \code{list("PCSpearman.adj",adjustment method)} has to be provided in the \code{methodlist} argument, where \code{adjustment method} has to be one of the options for multple testing adjustment provided by the standard \code{p.adjust} R function, i.e. one of \code{"BH"}, \code{"bonferroni"}, \code{"BY"}, \code{"fdr"}, \code{"hochberg"}, \code{"holm"} or \code{"hommel"}.}
#' \item{\code{list("DistCorr")} \cr Edge weights are estimated using distance correlation, where unadjusted p-values are employed to determine significance. To apply this method, the expression \code{list("DistCorr")} has to be provided in the \code{methodlist} argument. Note that the calculations may require larger computation times, as a permutation test is involved to derive the corresponding p-values for the distance correlations.}
#' \item{\code{list("DistCorr.adj",adjustment method)} \cr Edge weights are estimated using distance correlation, where p-values adjusted for multiple testing are employed to determine significance. To apply this method, the expression \code{list("DistCorr.adj",adjustment method)} has to be provided in the \code{methodlist} argument, where \code{adjustment method} has to be one of the options for multple testing adjustment provided by the standard \code{p.adjust} R function, i.e. one of \code{"BH"}, \code{"bonferroni"}, \code{"BY"}, \code{"fdr"}, \code{"hochberg"}, \code{"holm"} or \code{"hommel"}. Note that the calculations may require larger computation times, as a permutation test is involved to derive the corresponding p-values for the distance correlations.}
#' \item{\code{list("EBICglasso",correlation type,tuning parameter)} \cr Edge weights are estimated using the EBICglasso approach. To apply this method, the expression \code{list("EBICglasso",correlation type,tuning parameter)} has to be provided in the \code{methodlist} argument. Here, \code{correlation type} has to be one of the association options provided by the standard \code{cor} R function, i.e. one of \code{"kendall"}, \code{"pearson"} or \code{"spearman"}. Moreover, \code{tuning parameter} has to be a number specifying the EBIC tuning parameter \eqn{\gamma}. Typical choices include values between 0 and 0.5, where smaller values usually lead to a higher sensitivity in that more edges are included into the network. \cr Note that for EBICglasso, an additional specification of the \code{thresh} argument is obsolete, as it is not used for the application of the method.}
#' }
#' @return the adjacency matrix generated from the input data table \code{A} according to the specified estimation method
#' @examples
#' create.adjacency.matrix(ExDataA,methodlist=list("Spearman"))
#' create.adjacency.matrix(ExDataA,methodlist=list("Spearman.adj","bonferroni"))
#' create.adjacency.matrix(ExDataA,methodlist=list("PCSpearman"),thresh=0.1)
#' create.adjacency.matrix(ExDataA,methodlist=list("PCSpearman.adj","BH"),thresh=0.1)
#' create.adjacency.matrix(ExDataA,methodlist=list("DistCorr"))
#' create.adjacency.matrix(ExDataA,methodlist=list("DistCorr.adj","bonferroni"))
#' create.adjacency.matrix(ExDataB,methodlist=list("EBICglasso","spearman",0.1))
#' create.adjacency.matrix(ExDataB,methodlist=list("EBICglasso","pearson",0.05))
#' @export
#'
create.adjacency.matrix<-function (A, methodlist, thresh = 0.05)
  {
  stopifnot(`A needs to be a data frame, array, tibble or matrix` = (class(A) %in% c("tbl_df", "tbl", "data.frame","matrix", "array")), `The matrix A needs to have more than 4 rows` = nrow(A) >
              4, `The matrix A needs to have at least one column` = ncol(A) >=
              1, `Methodlist needs to be a list of Strings  with at least one element (except for the method EBICglasso where the third elements needs to be a number)` = class(methodlist) ==
              "list" && length(methodlist) >= 1 && (length(methodlist) ==
                                                      1 || class(methodlist[[2]]) == "character") ||
              methodlist[[1]] == "EBICglasso" && class(methodlist[[2]]) ==
              "character" && (class(methodlist[[3]]) == "numeric" ||
                                class(methodlist[[3]]) == "integer"), `For methods with the ending '.adj' a second method is needed and when using the method 'EBICglasso' one extra method and one extra value are needed` = (stringr::str_detect(methodlist[[1]],
                                                                                                                                                                                                                                                   "\\.adj") == FALSE && length(methodlist) >= 1) &&
              methodlist[[1]] != "EBICglasso" || (stringr::str_detect(methodlist[[1]],
                                                                      "\\.adj") && length(methodlist) >= 2) || (methodlist[[1]] ==
                                                                                                                  "EBICglasso" && length(methodlist) >= 3), `Thresh needs to be a non-negative number` = class(thresh) ==
              "numeric" && thresh >= 0)
    method <- methodlist[[1]]
    if (method == "Spearman") {
      correls <- Hmisc::rcorr(as.matrix(A), type = "spearman")
      cm <- correls$r
      diag(cm) <- 0
      cm[correls$P > thresh] <- 0
    }
    if (method == "PCSpearman") {
      correls <- ppcor::pcor(as.matrix(A), method = "spearman")
      cm <- correls$estimate
      diag(cm) <- 0
      cm[correls$p.value > thresh] <- 0
    }
    if (method == "Spearman.adj") {
      correls <- Hmisc::rcorr(as.matrix(A), type = "spearman")
      cm <- correls$r
      pmat <- correls$P
      ltri <- lower.tri(pmat)
      utri <- upper.tri(pmat)
      adj.method <- methodlist[[2]]
      pmat[ltri] <- stats::p.adjust(pmat[ltri], method = adj.method)
      pmat[utri] <- t(pmat)[utri]
      diag(cm) <- 0
      cm[pmat > thresh] <- 0
    }
    if (method == "PCSpearman.adj") {
      correls <- ppcor::pcor(as.matrix(A), method = "spearman")
      cm <- correls$estimate
      pmat <- correls$p.value
      ltri <- lower.tri(pmat)
      utri <- upper.tri(pmat)
      adj.method <- methodlist[[2]]
      pmat[ltri] <- stats::p.adjust(pmat[ltri], method = adj.method)
      pmat[utri] <- t(pmat)[utri]
      diag(cm) <- 0
      cm[pmat > thresh] <- 0
    }
    if (method == "DistCorr") {
      A <- as.matrix(A)
      cm <- array(data = NA, dim = c(dim(A)[2], dim(A)[2]))
      for (i in 1:dim(A)[2]) {
        for (j in 1:dim(A)[2]) {
          cm[i, j] <- energy::dcor(A[, i], A[, j], index = 1)
        }
      }
      rownames(cm) <- colnames(A)
      colnames(cm) <- colnames(A)
      pvals <- array(data = NA, dim = c(dim(A)[2], dim(A)[2]))
      for (i in 1:dim(A)[2]) {
        for (j in 1:dim(A)[2]) {
          set.seed(24)
          pvals[i, j] <- energy::dcor.test(A[, i], A[, j], index = 1,
                                           R = 10000)$p.value
        }
      }
      rownames(pvals) <- colnames(A)
      colnames(pvals) <- colnames(A)
      diag(cm) <- 0
      cm[pvals > thresh] <- 0
    }
    if (method == "DistCorr.adj") {
      A <- as.matrix(A)
      cm <- array(data = NA, dim = c(dim(A)[2], dim(A)[2]))
      for (i in 1:dim(A)[2]) {
        for (j in 1:dim(A)[2]) {
          cm[i, j] <- energy::dcor(A[, i], A[, j], index = 1)
        }
      }
      rownames(cm) <- colnames(A)
      colnames(cm) <- colnames(A)
      pvals <- array(data = NA, dim = c(dim(A)[2], dim(A)[2]))
      for (i in 1:dim(A)[2]) {
        for (j in 1:dim(A)[2]) {
          set.seed(24)
          pvals[i, j] <- energy::dcor.test(A[, i], A[, j], index = 1,
                                           R = 10000)$p.value
        }
      }
      rownames(pvals) <- colnames(A)
      colnames(pvals) <- colnames(A)
      pmat <- pvals
      ltri <- lower.tri(pmat)
      utri <- upper.tri(pmat)
      adj.method <- methodlist[[2]]
      pmat[ltri] <- stats::p.adjust(pmat[ltri], method = adj.method)
      pmat[utri] <- t(pmat)[utri]
      rownames(pmat) <- colnames(A)
      colnames(pmat) <- colnames(A)
      diag(cm) <- 0
      cm[pmat > thresh] <- 0
    }
    if (method == "EBICglasso") {
      corm.ebicglasso <- methodlist[[2]]
      tun.ebicglasso <- methodlist[[3]]
      correls <- stats::cor(A, method = corm.ebicglasso)
      cm <- qgraph::EBICglasso(correls, n = dim(A)[1], gamma = tun.ebicglasso)
      diag(cm) <- 0
    }
    output <- cm
    return(output)
  }


#' Creates \code{igraph} graph and MST results
#' @description
#' This function creates \code{igraph} graph and corresponding minimum spanning tree (MST) results from an adjacency matrix that is produced out of an input data table, using an estimation method specified by the user.
#' @param A input data table from which the adjacency matrix will be generated, to be provided in form of a matrix, array, data frame or tibble
#' @param methodlist a list specifying the method which is used to estimate and create the adjacency matrix; see details for further information
#' @param thresh a number between 0 and 1 (default is set to 0.05) specifying the singificance level: if the p-value corresponding to an edge weight is greater than \code{thresh}, the corresponding edge weight is not considered to be significant and thus set to zero
#' @details
#' This function creates \code{igraph} graph and corresponding minimum spanning tree (MST; derived using Prim's algorithm) results from an adjacency matrix that is produced out of an input data table, using an estimation method specified by the user. The function builds on respective implementations in the \code{igraph} R package.
#' \cr
#' The network estimation method has to be specified in form of a list in the \code{methodlist} argument. Currently, the following estimation methods are supported:
#' \itemize{
#' \item{\code{list("Spearman")} \cr Edge weights are estimated using Spearman correlation, where unadjusted p-values are employed to determine significance. To apply this method, the expression \code{list("Spearman")} has to be provided in the \code{methodlist} argument.}
#' \item{\code{list("Spearman.adj",adjustment method)} \cr Edge weights are estimated using Spearman correlation, where p-values adjusted for multiple testing are employed to determine significance. To apply this method, the expression \code{list("Spearman.adj",adjustment method)} has to be provided in the \code{methodlist} argument, where \code{adjustment method} has to be one of the options for multple testing adjustment provided by the standard \code{p.adjust} R function, i.e. one of \code{"BH"}, \code{"bonferroni"}, \code{"BY"}, \code{"fdr"}, \code{"hochberg"}, \code{"holm"} or \code{"hommel"}.}
#' \item{\code{list("PCSpearman")} \cr Edge weights are estimated using partial Spearman correlation, where unadjusted p-values are employed to determine significance. To apply this method, the expression \code{list("PCSpearman")} has to be provided in the \code{methodlist} argument.}
#' \item{\code{list("PCSpearman.adj",adjustment method)} \cr Edge weights are estimated using partial Spearman correlation, where p-values adjusted for multiple testing are employed to determine significance. To apply this method, the expression \code{list("PCSpearman.adj",adjustment method)} has to be provided in the \code{methodlist} argument, where \code{adjustment method} has to be one of the options for multple testing adjustment provided by the standard \code{p.adjust} R function, i.e. one of \code{"BH"}, \code{"bonferroni"}, \code{"BY"}, \code{"fdr"}, \code{"hochberg"}, \code{"holm"} or \code{"hommel"}.}
#' \item{\code{list("DistCorr")} \cr Edge weights are estimated using distance correlation, where unadjusted p-values are employed to determine significance. To apply this method, the expression \code{list("DistCorr")} has to be provided in the \code{methodlist} argument. Note that the calculations may require larger computation times, as a permutation test is involved to derive the corresponding p-values for the distance correlations.}
#' \item{\code{list("DistCorr.adj",adjustment method)} \cr Edge weights are estimated using distance correlation, where p-values adjusted for multiple testing are employed to determine significance. To apply this method, the expression \code{list("DistCorr.adj",adjustment method)} has to be provided in the \code{methodlist} argument, where \code{adjustment method} has to be one of the options for multple testing adjustment provided by the standard \code{p.adjust} R function, i.e. one of \code{"BH"}, \code{"bonferroni"}, \code{"BY"}, \code{"fdr"}, \code{"hochberg"}, \code{"holm"} or \code{"hommel"}. Note that the calculations may require larger computation times, as a permutation test is involved to derive the corresponding p-values for the distance correlations.}
#' \item{\code{list("EBICglasso",correlation type,tuning parameter)} \cr Edge weights are estimated using the EBICglasso approach. To apply this method, the expression \code{list("EBICglasso",correlation type,tuning parameter)} has to be provided in the \code{methodlist} argument. Here, \code{correlation type} has to be one of the association options provided by the standard \code{cor} R function, i.e. one of \code{"kendall"}, \code{"pearson"} or \code{"spearman"}. Moreover, \code{tuning parameter} has to be a number specifying the EBIC tuning parameter \eqn{\gamma}. Typical choices include values between 0 and 0.5, where smaller values usually lead to a higher sensitivity in that more edges are included into the network. \cr Note that for EBICglasso, an additional specification of the \code{thresh} argument is obsolete, as it is not used for the application of the method.}
#' }
#' @return a list with 15 elements: the adjacency matrix, the \code{igraph} graph, the graph communities (derived by the Girvan-Newman algorithm based on edge betweenness), the graph clustering, the graph vertex degrees, the graph number of edges, the graph number of clusters, the graph number of isolated nodes, the minimum spanning tree (MST), the MST communities (derived by the Girvan-Newman algorithm based on edge betweenness), the MST clustering, the MST vertex degrees, the MST number of edges, the MST number of clusters and the MST number of isolated nodes
#' @examples
#' create.graph(ExDataA,methodlist=list("Spearman"))
#' create.graph(ExDataA,methodlist=list("Spearman.adj","bonferroni"))
#' create.graph(ExDataA,methodlist=list("PCSpearman"),thresh=0.1)
#' create.graph(ExDataA,methodlist=list("PCSpearman.adj","BH"),thresh=0.1)
#' create.graph(ExDataA,methodlist=list("DistCorr"))
#' create.graph(ExDataA,methodlist=list("DistCorr.adj","bonferroni"))
#' create.graph(ExDataB,methodlist=list("EBICglasso","spearman",0.1))
#' create.graph(ExDataB,methodlist=list("EBICglasso","pearson",0.05))
#' @export
#'
create.graph<-function (A, methodlist, thresh = 0.05)
  {
    cm <- create.adjacency.matrix(A, methodlist, thresh)
    if (all(cm == 0)) {
      output <- list(cm, NA, NA, NA, NA, 0, ncol(cm), ncol(cm), NA, NA,
                     NA, NA, 0, ncol(cm), ncol(cm))
      names(output) <- c("adjacency matrix", "graph",
                         "graph communities", "graph clustering",
                         "graph vertex degrees", "graph number of edges",
                         "graph number of clusters", "graph number of isolated nodes",
                         "MST", "MST communities", "MST clustering",
                         "MST vertex degrees", "MST number of edges",
                         "MST number of clusters", "MST number of isolated nodes")
      return(output)
    }
    else {
      g <- igraph::graph.adjacency(cm, weighted = TRUE, mode = "undirected",
                                   diag = FALSE)
      g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
      igraph::E(g)[which(igraph::E(g)$weight < 0)]$color <- "red"
      igraph::E(g)[which(igraph::E(g)$weight > 0)]$color <- "blue"
      igraph::V(g)$name <- igraph::V(g)$name
      g.communities <- igraph::edge.betweenness.community(g, weights = NA,
                                                          directed = FALSE)
      g.clustering <- igraph::make_clusters(g, membership = g.communities$membership)
      igraph::V(g)$color <- g.communities$membership
      deg.g <- igraph::degree(g)
      edgenum.g <- igraph::gsize(g)
      clusnum.g <- length(unique(g.communities$membership))
      numisolnodes.g <- length(which(deg.g == 0))
      mst <- igraph::mst(g, algorithm = "prim")
      mst.communities <- igraph::edge.betweenness.community(mst, weights = NA,
                                                            directed = FALSE)
      mst.clustering <- igraph::make_clusters(mst, membership = mst.communities$membership)
      igraph::V(mst)$color <- mst.communities$membership
      deg.mst <- igraph::degree(mst)
      edgenum.mst <- igraph::gsize(mst)
      clusnum.mst <- length(unique(mst.communities$membership))
      numisolnodes.mst <- length(which(deg.mst == 0))
      output <- list(cm, g, g.communities, g.clustering, deg.g,
                     edgenum.g, clusnum.g, numisolnodes.g, mst, mst.communities,
                     mst.clustering, deg.mst, edgenum.mst, clusnum.mst,
                     numisolnodes.mst)
      names(output) <- c("adjacency matrix", "graph",
                         "graph communities", "graph clustering",
                         "graph vertex degrees", "graph number of edges",
                         "graph number of clusters", "Graph number of isolated nodes",
                         "MST", "MST communities", "MST clustering",
                         "MST vertex degrees", "MST number of edges",
                         "MST number of clusters", "MST number of isolated nodes")
      return(output)
    }
  }



###2.) Network differences

######(i) overall network difference characteristics

#' Calculates the Frobenius metric
#' @description This function calculates the Frobenius metric between two adjacency matrices of the same dimension.
#' @param A,B adjacency matrices of the same dimension
#' @details
#' This function calculates the Frobenius metric between two adjacency matrices of the same dimension. The Frobenius metric is an overall characteristic that can be employed to compare two networks.
#' @return the Frobenius metric between the two adjacency matrices \code{A} and \code{B}
#' @examples
#' A<-create.adjacency.matrix(ExDataA,methodlist=list("Spearman"))
#' B<-create.adjacency.matrix(ExDataB,methodlist=list("Spearman"))
#' frobenius.metric(A,B)
#' @export
#'
frobenius.metric<-function(A,B){
  stopifnot(`A and B need to be a data frame, array, tibble or matrix` = all((class(A) %in% c("tbl_df", "tbl", "data.frame","matrix", "array")) &
                                                                         (class(B) %in% c("tbl_df", "tbl", "data.frame","matrix", "array"))),
            `A and B need to have the same dimensions` = all(nrow(A)==nrow(B) & ncol(A)==ncol(B)))

  output<-sqrt(sum(abs((A-B)^2)))
  return(output)
}

#' Calculates the maximum metric
#' @description This function calculates the maximum metric between two adjacency matrices of the same dimension.
#' @param A,B adjacency matrices of the same dimension
#' @details
#' This function calculates the maximum metric between two adjacency matrices of the same dimension. The maximum metric is an overall characteristic that can be employed to compare two networks.
#' @return the maximum metric between the two adjacency matrices \code{A} and \code{B}
#' @examples
#' A<-create.adjacency.matrix(ExDataA,methodlist=list("Spearman"))
#' B<-create.adjacency.matrix(ExDataB,methodlist=list("Spearman"))
#' maximum.metric(A,B)
#' @export
#'
maximum.metric<-function(A,B){
  stopifnot(`A and B need to be a data frame, array, tibble or matrix` = all((class(A) %in% c("tbl_df", "tbl", "data.frame","matrix", "array")) &
                                                                         (class(B) %in% c("tbl_df", "tbl", "data.frame","matrix", "array"))),
            `A and B need to have the same dimensions` = all(nrow(A)==nrow(B) & ncol(A)==ncol(B)))

  output<-max(abs(A-B))
  return(output)
}


#' Calculates the spectral distance
#' @description This function calculates the spectral distance between two adjacency matrices of the same dimension.
#' @param A,B adjacency matrices of the same dimension
#' @details
#' This function calculates the spectral distance between two adjacency matrices of the same dimension. The spectral distance is an overall characteristic that can be employed to compare two networks.
#' @return the spectral distance between the two adjacency matrices \code{A} and \code{B}
#' @examples
#' A<-create.adjacency.matrix(ExDataA,methodlist=list("Spearman"))
#' B<-create.adjacency.matrix(ExDataB,methodlist=list("Spearman"))
#' spec.dist(A,B)
#' @export
#'
spec.dist<-function(A,B){
  stopifnot(`A and B need to be a data frame, array, tibble or matrix` = all((class(A) %in% c("tbl_df", "tbl", "data.frame","matrix", "array")) &
                                                                         (class(B) %in% c("tbl_df", "tbl", "data.frame","matrix", "array"))),
            `A and B need to have the same dimensions` = all(nrow(A)==nrow(B) & ncol(A)==ncol(B)))

  output<-sqrt(sum(((eigen(A)$values)-(eigen(B)$values))^2))
  return(output)
}


#' Calculates the Jaccard distance
#' @description This function calculates the Jaccard distance between two adjacency matrices of the same dimension.
#' @param A,B adjacency matrices of the same dimension
#' @details
#' This function calculates the Jaccard distance between two adjacency matrices of the same dimension. The Jaccard distance is an overall characteristic that can be employed to compare two networks.
#' @return the Jaccard distance between the two adjacency matrices \code{A} and \code{B}
#' @examples
#' A<-create.adjacency.matrix(ExDataA,methodlist=list("Spearman"))
#' B<-create.adjacency.matrix(ExDataB,methodlist=list("Spearman"))
#' jaccard.dist(A,B)
#' @export
#'
jaccard.dist<-function(A,B){
  stopifnot(`A and B need to be a data frame, array, tibble or matrix` = all((class(A) %in% c("tbl_df", "tbl", "data.frame","matrix", "array")) &
              (class(B) %in% c("tbl_df", "tbl", "data.frame","matrix", "array"))),
            `A and B need to have the same dimensions` = all(nrow(A)==nrow(B) & ncol(A)==ncol(B)))


  A.new <- data.frame(var1 = rownames(A)[row(A)[upper.tri(A)]],
                      var2 = colnames(A)[col(A)[upper.tri(A)]],
                      corr.val = A[upper.tri(A)])
  idx.A<-which(A.new[,"corr.val"]!=0)
  edges.A<-paste0(paste0(A.new[idx.A,"var1"],"-"),A.new[idx.A,"var2"])


  B.new <- data.frame(var1 = rownames(B)[row(B)[upper.tri(B)]],
                      var2 = colnames(B)[col(B)[upper.tri(B)]],
                      corr.val = B[upper.tri(B)])
  idx.B<-which(B.new[,"corr.val"]!=0)
  edges.B<-paste0(paste0(B.new[idx.B,"var1"],"-"),B.new[idx.B,"var2"])

  len.intersect.AB<-length(intersect(edges.A,edges.B))
  len.union.AB<-length(edges.A)+length(edges.B)-len.intersect.AB
  jaccard.index<-len.intersect.AB/len.union.AB

  output<-1-jaccard.index
  return(output)
}



#' Calculates the difference with respect to global strength
#' @description This function calculates the difference with respect to the global strength between two adjacency matrices of the same dimension.
#' @param A,B adjacency matrices of the same dimension
#' @details
#' This function calculates the difference with respect to the global strength between two adjacency matrices of the same dimension. The global strength is an overall characteristic that can be employed to compare two networks.
#' @return the difference with respect to the global strength between the two adjacency matrices \code{A} and \code{B}.
#' @examples
#' A<-create.adjacency.matrix(ExDataA,methodlist=list("Spearman"))
#' B<-create.adjacency.matrix(ExDataB,methodlist=list("Spearman"))
#' global.str(A,B)
#' @export
#'
global.str<-function(A,B){
  stopifnot(`A and B need to be a data frame, array, tibble or matrix` = all((class(A) %in% c("tbl_df", "tbl", "data.frame","matrix", "array")) &
                                                                         (class(B) %in% c("tbl_df", "tbl", "data.frame","matrix", "array"))),
            `A and B need to have the same dimensions` = all(nrow(A)==nrow(B) & ncol(A)==ncol(B)))

  output<-abs(sum(abs(A))-sum(abs(B)))
  return(output)
}


#' Calculates the differences with respect to the number of edges, clusters and isolated nodes between two networks
#' @description This function calculates the differences with respect to the number of edges, clusters (obtained by the Girvan-Newman algorithm) and isolated nodes (i.e. nodes without any edge) between two networks, for both \code{igraph} graphs and the corresponding minimum spanning trees (MSTs).
#' @param A,B output results when applying the function \code{create.graph} to two data tables from which the two networks are constructed
#' @details
#' This function calculates the differences with respect to the number of edges, clusters (obtained by the Girvan-Newman algorithm) and isolated nodes (i.e. nodes without any edge) between two networks, for both \command{igraph} graphs and the corresponding minimum spanning trees (MSTs). It builds on the \command{create.graph} function, which creates \command{igraph} graph and corresponding MST results from adjacency matrices that are produced out of input tables, using an estimation method specified by the user; see the documentation of \command{create.graph} for further information. Differences in numbers of edges, clusters and isolated nodes are overall characteristics that can be employed to compare two networks.
#' @return a vector of length 6 containing the graph difference in the number of edges, the graph difference in the number of clusters, the graph difference in the number of isolated nodes, the MST difference in the number of edges, the MST difference in the number of clusters and the MST difference in the number of isolated nodes
#' @examples
#' A<-create.graph(ExDataA,methodlist=list("Spearman"))
#' B<-create.graph(ExDataB,methodlist=list("Spearman"))
#' number.differences(A,B)
#' @export
#'
number.differences<-function(A,B){

  stopifnot(`A and B need to be adjacency matrices.` = all(class(A[[1]]) %in% c("matrix","array") & class(B[[1]]) %in% c("matrix","array")
            & nrow(A[[1]])==ncol(A[[1]]) & nrow(B[[1]])==ncol(B[[1]])),
            `The listelements 6, 7, 8, 13, 14 and 15 need to be numbers. (This is the case if A and B are results of the function create.graph.)` =
            all(all(lapply(A, class)[6:8] %in% c("numeric", "integer")) &
              all(lapply(B, class)[6:8] %in% c("numeric", "integer")) &
              all(lapply(A, class)[13:15] %in% c("numeric", "integer")) &
              all(lapply(B, class)[13:15] %in% c("numeric", "integer")))
  )


    dnumedges.g<-abs(A[[6]]-B[[6]])
    dnumclus.g<-abs(A[[7]]-B[[7]])
    dnumisolnodes.g<-abs(A[[8]]-B[[8]])


    dnumedges.mst<-abs(A[[13]]-B[[13]])
    dnumclus.mst<-abs(A[[14]]-B[[14]])
    dnumisolnodes.mst<-abs(A[[15]]-B[[15]])


    output<-c(dnumedges.g,dnumclus.g,dnumisolnodes.g,dnumedges.mst,dnumclus.mst,dnumisolnodes.mst)
    names(output)<-c("Graph Diff num of edges","Graph Diff num of clusters","Graph diff num of isolated nodes","MST Diff num of edges","MST Diff num of clusters","MST diff num of isolated nodes")
    return(output)
  }




######(ii) node-specifc network difference characteristics

#' Calculates differences in degree for each network node
#' @description This function calculates the differences in degree for each node between two networks (adjacency matrices of the same dimension).
#' @param X,Y adjacency matrices of the same dimension
#' @details
#' This function calculates the differences in degree for each node between two networks (adjacency matrices of the same dimension). Differences in degree are one of the node-specific network difference characteristics to compare two networks.
#' @return a vector of length \eqn{N} (number of nodes), containing the differences in degree between the two networks for each node
#' @examples
#' X<-create.adjacency.matrix(ExDataA,methodlist=list("Spearman"))
#' Y<-create.adjacency.matrix(ExDataB,methodlist=list("Spearman"))
#' degree.inv(X,Y)
#' @export
#'
degree.inv<-function(X,Y){

  stopifnot(`X and Y need to be matrices` = all(class(X) %in% c("matrix","array") & class(Y) %in% c("matrix","array")),
            `X and Y need to have the same dimensions` = all(nrow(X)==nrow(Y) & ncol(X)==ncol(Y)),
            `To compare the graphs of X and Y correctly they need the same columntitles` = colnames(X)==colnames(Y))


  if(all(X==0)){
    A<-igraph::graph.adjacency(X,weighted=TRUE, mode="undirected", diag=FALSE)
  }else{
    A<-igraph::graph.adjacency(X,weighted=TRUE, mode="undirected", diag=FALSE)
    # Simplfy the adjacency object
    A<-igraph::simplify(A, remove.multiple=TRUE, remove.loops=TRUE)
    # Convert edge weights to absolute values
    igraph::E(A)$weight <- abs(igraph::E(A)$weight)
  }

  if(all(Y==0)){
    B<-igraph::graph.adjacency(Y,weighted=TRUE, mode="undirected", diag=FALSE)
  }else{
    B<-igraph::graph.adjacency(Y,weighted=TRUE, mode="undirected", diag=FALSE)
    # Simplfy the adjacency object
    B<-igraph::simplify(B, remove.multiple=TRUE, remove.loops=TRUE)
    # Convert edge weights to absolute values
    igraph::E(B)$weight <- abs(igraph::E(B)$weight)
  }

  degA<-igraph::degree(A)
  degB<-igraph::degree(B)

  output<-abs(degA-degB)
  names(output)<-colnames(X)

  return(output)
}


#' Calculates differences in betweenness for each network node
#' @description This function calculates the differences in betweenness for each node between two networks (adjacency matrices of the same dimension).
#' @param X,Y adjacency matrices of the same dimension
#' @details
#' This function calculates the differences in betweenness for each node between two networks (adjacency matrices of the same dimension). Differences in betweenness are one of the node-specific network difference characteristics to compare two networks.
#' @return a vector of length \eqn{N} (number of nodes), containing the differences in betweenness between the two networks for each node
#' @examples
#' X<-create.adjacency.matrix(ExDataA,methodlist=list("Spearman"))
#' Y<-create.adjacency.matrix(ExDataB,methodlist=list("Spearman"))
#' betweenness.inv(X,Y)
#' @export
#'
betweenness.inv<-function(X,Y){

  stopifnot(`X and Y need to be matrices` = all(class(X) %in% c("matrix","array") & class(Y) %in% c("matrix","array")),
            `X and Y need to have the same dimensions` = all(nrow(X)==nrow(Y) & ncol(X)==ncol(Y)),
            `To compare the graphs of X and Y correctly they need the same columntitles` = colnames(X)==colnames(Y))

  if(all(X==0)){
    A<-igraph::graph.adjacency(X,weighted=TRUE, mode="undirected", diag=FALSE)
  }else{
    A<-igraph::graph.adjacency(X,weighted=TRUE, mode="undirected", diag=FALSE)
    # Simplfy the adjacency object
    A<-igraph::simplify(A, remove.multiple=TRUE, remove.loops=TRUE)
    # Convert edge weights to absolute values
    igraph::E(A)$weight <- abs(igraph::E(A)$weight)
  }

  if(all(Y==0)){
    B<-igraph::graph.adjacency(Y,weighted=TRUE, mode="undirected", diag=FALSE)
  }else{
    B<-igraph::graph.adjacency(Y,weighted=TRUE, mode="undirected", diag=FALSE)
    # Simplfy the adjacency object
    B<-igraph::simplify(B, remove.multiple=TRUE, remove.loops=TRUE)
    # Convert edge weights to absolute values
    igraph::E(B)$weight <- abs(igraph::E(B)$weight)
  }

  betA<-igraph::betweenness(A)
  betB<-igraph::betweenness(B)

  output<-abs(betA-betB)
  names(output)<-colnames(X)

  return(output)
}


#' Calculates differences in closeness for each network node
#' @description This function calculates the differences in closeness for each node between two networks (adjacency matrices of the same dimension).
#' @param X,Y adjacency matrices of the same dimension
#' @details
#' This function calculates the differences in closeness for each node between two networks (adjacency matrices of the same dimension). Differences in closeness are one of the node-specific network difference characteristics to compare two networks.
#' @return a vector of length \eqn{N} (number of nodes), containing the differences in closeness between the two networks for each node
#' @examples
#' X<-create.adjacency.matrix(ExDataA,methodlist=list("Spearman"))
#' Y<-create.adjacency.matrix(ExDataB,methodlist=list("Spearman"))
#' closeness.inv(X,Y)
#' @export
#'
closeness.inv<-function(X,Y){

  stopifnot(`X and Y need to be matrices` = all(class(X) %in% c("matrix","array") & class(Y) %in% c("matrix","array")),
            `X and Y need to have the same dimensions` = all(nrow(X)==nrow(Y) && ncol(X)==ncol(Y)),
            `To compare the graphs of X and Y correctly they need the same columntitles` = colnames(X)==colnames(Y))


  if(all(X==0)){
    A<-igraph::graph.adjacency(X,weighted=TRUE, mode="undirected", diag=FALSE)
  }else{
    A<-igraph::graph.adjacency(X,weighted=TRUE, mode="undirected", diag=FALSE)
    # Simplify the adjacency object
    A<-igraph::simplify(A, remove.multiple=TRUE, remove.loops=TRUE)
    # Convert edge weights to absolute values
    igraph::E(A)$weight <- abs(igraph::E(A)$weight)
  }

  if(all(Y==0)){
    B<-igraph::graph.adjacency(Y,weighted=TRUE, mode="undirected", diag=FALSE)
  }else{
    B<-igraph::graph.adjacency(Y,weighted=TRUE, mode="undirected", diag=FALSE)
    # Simplfy the adjacency object
    B<-igraph::simplify(B, remove.multiple=TRUE, remove.loops=TRUE)
    # Convert edge weights to absolute values
    igraph::E(B)$weight <- abs(igraph::E(B)$weight)
  }

  cloA<-igraph::closeness(A)
  cloB<-igraph::closeness(B)

  output<-abs(cloA-cloB)
  names(output)<-colnames(X)

  return(output)
}


#' Calculates differences in eigenvector centrality for each network node
#' @description This function calculates the differences in eigenvector centrality for each node between two networks (adjacency matrices of the same dimension).
#' @param X,Y adjacency matrices of the same dimension
#' @details
#' This function calculates the differences in eigenvector centrality for each node between two networks (adjacency matrices of the same dimension). Differences in eigenvector centrality are one of the node-specific network difference characteristics to compare two networks.
#' @return a vector of length \eqn{N} (number of nodes), containing the differences in eigenvector centrality between the two networks for each node
#' @examples
#' X<-create.adjacency.matrix(ExDataA,methodlist=list("Spearman"))
#' Y<-create.adjacency.matrix(ExDataB,methodlist=list("Spearman"))
#' eigen.inv(X,Y)
#' @export
#'
eigen.inv<-function(X,Y){

  stopifnot(`X and Y need to be matrices` = all(class(X) %in% c("matrix","array") & class(Y) %in% c("matrix","array")),
            `X and Y need to have the same dimensions` = all(nrow(X)==nrow(Y) && ncol(X)==ncol(Y)),
            `To compare the graphs of X and Y correctly they need the same columntitles` = colnames(X)==colnames(Y))

  if(all(X==0)){
    A<-igraph::graph.adjacency(X,weighted=TRUE, mode="undirected", diag=FALSE)
  }else{
    A<-igraph::graph.adjacency(X,weighted=TRUE, mode="undirected", diag=FALSE)
    # Simplfy the adjacency object
    A<-igraph::simplify(A, remove.multiple=TRUE, remove.loops=TRUE)
    # Convert edge weights to absolute values
    igraph::E(A)$weight <- abs(igraph::E(A)$weight)
  }

  if(all(Y==0)){
    B<-igraph::graph.adjacency(Y,weighted=TRUE, mode="undirected", diag=FALSE)
  }else{
    B<-igraph::graph.adjacency(Y,weighted=TRUE, mode="undirected", diag=FALSE)
    # Simplfy the adjacency object
    B<-igraph::simplify(B, remove.multiple=TRUE, remove.loops=TRUE)
    # Convert edge weights to absolute values
    igraph::E(B)$weight <- abs(igraph::E(B)$weight)
  }

  ecA<-igraph::eigen_centrality(A)$vector
  ecB<-igraph::eigen_centrality(B)$vector

  output<-abs(ecA-ecB)
  names(output)<-colnames(X)

  return(output)
}


######(iii) edge-specific network difference characteristics

##a. taking account of directions (signs) of correlations

#' Calculates the differences in edge weight for each network edge
#' @description This function calculates the differences in edge weight for each edge between two networks (adjacency matrices of the same dimension).
#' @param A,B adjacency matrices of the same dimension
#' @details
#'This function calculates the differences in edge weight for each edge between two networks (adjacency matrices of the same dimension). In particular, it takes account of potentially different directions (signs) of the edge weights (associations) when deriving the differences. Differences in edge weight are one of the edge-specific network difference characteristics to compare two networks.
#' @return a symmetric \eqn{N \times N} matrix, with \eqn{N} denoting the number of nodes, containing the differences in edge weight between two respective nodes
#' @examples
#' A<-create.adjacency.matrix(ExDataA,methodlist=list("Spearman"))
#' B<-create.adjacency.matrix(ExDataB,methodlist=list("Spearman"))
#' edge.inv.direc(A,B)
#' @export
#'
edge.inv.direc<-function(A,B){

  stopifnot(`A and B need to be matrices or tibbles` = all(class(A) %in% c("matrix","array","tbl_df","tbl","data.frame") &
                                                       class(B) %in% c("matrix","array","tbl_df","tbl","data.frame")),
            `A and B need to have the same dimensions` = all(nrow(A)==nrow(B) && ncol(A)==ncol(B)))

  output<-abs(A-B)
  rownames(output)<-rownames(A)
  colnames(output)<-colnames(A)
  return(output)
}


##b. not taking account of directions (signs) of correlations

#'Calculates the differences in absolute edge weight for each network edge
#' @description This function calculates the differences in absolute edge weight for each edge between two networks (adjacency matrices of the same dimension).
#' @param A,B adjacency matrices of the same dimension
#' @details
#' This function calculates the differences in absolute edge weight for each edge between two networks (adjacency matrices of the same dimension). In particular, by considering absolute edge weights, it does not take account of directions (signs) of the edge weights (associations) when deriving the differences. Differences in absolute edge weight are one of the edge-specific network difference characteristics to compare two networks.
#' @return a symmetric \eqn{N \times N} matrix, with \eqn{N} denoting the number of nodes, containing the differences in absolute edge weight between two respective nodes
#' @examples
#' A<-create.adjacency.matrix(ExDataA,methodlist=list("Spearman"))
#' B<-create.adjacency.matrix(ExDataB,methodlist=list("Spearman"))
#' edge.inv(A,B)
#' @export
#'
edge.inv<-function(A,B){

  stopifnot(`A and B need to be matrices or tibbles` = all(class(A) %in% c("matrix","array","tbl_df","tbl","data.frame") &
                                                       class(B) %in% c("matrix","array","tbl_df","tbl","data.frame")),
            `A and B need to have the same dimensions` = all(nrow(A)==nrow(B) && ncol(A)==ncol(B)))

  output<-abs(abs(A)-abs(B))
  rownames(output)<-rownames(A)
  colnames(output)<-colnames(A)
  return(output)
}

###3.) Permutation test

#' Permutation-based test for differences between two networks
#' @description This function provides a permutation-based frame for testing for differences between two networks. In particular, various (i) network estimation methods and (ii) network difference characteristics can be specified.
#' @param A,B input data tables from which the adjacency matrices will be generated, to be provided in form of matrices, arrays, data frames or tibbles; need to have the same number of columns (corresponding to the number of nodes)
#' @param permnum a number, specifying the number of permutations
#' @param methodlist a list specifying the method which is used to estimate and create the adjacency matrices; see details for possible options and further information
#' @param thresh a number between 0 and 1 (default is set to 0.05) specifying the singificance level: if the p-value corresponding to an edge weight is greater than \code{thresh}, the corresponding edge weight is not considered to be significant and thus set to zero
#' @param score.funct the function used to compare the adjacency matrices A and B; see details for possible options and further information
#' @param paired Boolean, specifying whether the data underlying the two networks is paired or not
#' @details
#' This function provides a permutation-based frame for testing for differences between two networks. In particular, various (i) network estimation methods and (ii) network difference characteristics can be specified.
#' \cr
#' (i) The network estimation method has to be specified in form of a list in the \code{methodlist} argument. Currently, the following estimation methods are supported:
#' \itemize{
#' \item{\code{list("Spearman")} \cr Edge weights are estimated using Spearman correlation, where unadjusted p-values are employed to determine significance. To apply this method, the expression \code{list("Spearman")} has to be provided in the \code{methodlist} argument.}
#' \item{\code{list("Spearman.adj",adjustment method)} \cr Edge weights are estimated using Spearman correlation, where p-values adjusted for multiple testing are employed to determine significance. To apply this method, the expression \code{list("Spearman.adj",adjustment method)} has to be provided in the \code{methodlist} argument, where \code{adjustment method} has to be one of the options for multple testing adjustment provided by the standard \code{p.adjust} R function, i.e. one of \code{"BH"}, \code{"bonferroni"}, \code{"BY"}, \code{"fdr"}, \code{"hochberg"}, \code{"holm"} or \code{"hommel"}.}
#' \item{\code{list("PCSpearman")} \cr Edge weights are estimated using partial Spearman correlation, where unadjusted p-values are employed to determine significance. To apply this method, the expression \code{list("PCSpearman")} has to be provided in the \code{methodlist} argument.}
#' \item{\code{list("PCSpearman.adj",adjustment method)} \cr Edge weights are estimated using partial Spearman correlation, where p-values adjusted for multiple testing are employed to determine significance. To apply this method, the expression \code{list("PCSpearman.adj",adjustment method)} has to be provided in the \code{methodlist} argument, where \code{adjustment method} has to be one of the options for multple testing adjustment provided by the standard \code{p.adjust} R function, i.e. one of \code{"BH"}, \code{"bonferroni"}, \code{"BY"}, \code{"fdr"}, \code{"hochberg"}, \code{"holm"} or \code{"hommel"}.}
#' \item{\code{list("DistCorr")} \cr Edge weights are estimated using distance correlation, where unadjusted p-values are employed to determine significance. To apply this method, the expression \code{list("DistCorr")} has to be provided in the \code{methodlist} argument. Note that the calculations may require larger computation times, as a permutation test is involved to derive the corresponding p-values for the distance correlations.}
#' \item{\code{list("DistCorr.adj",adjustment method)} \cr Edge weights are estimated using distance correlation, where p-values adjusted for multiple testing are employed to determine significance. To apply this method, the expression \code{list("DistCorr.adj",adjustment method)} has to be provided in the \code{methodlist} argument, where \code{adjustment method} has to be one of the options for multple testing adjustment provided by the standard \code{p.adjust} R function, i.e. one of \code{"BH"}, \code{"bonferroni"}, \code{"BY"}, \code{"fdr"}, \code{"hochberg"}, \code{"holm"} or \code{"hommel"}. Note that the calculations may require larger computation times, as a permutation test is involved to derive the corresponding p-values for the distance correlations.}
#' \item{\code{list("EBICglasso",correlation type,tuning parameter)} \cr Edge weights are estimated using the EBICglasso approach. To apply this method, the expression \code{list("EBICglasso",correlation type,tuning parameter)} has to be provided in the \code{methodlist} argument. Here, \code{correlation type} has to be one of the association options provided by the standard \code{cor} R function, i.e. one of \code{"kendall"}, \code{"pearson"} or \code{"spearman"}. Moreover, \code{tuning parameter} has to be a number specifying the EBIC tuning parameter \eqn{\gamma}. Typical choices include values between 0 and 0.5, where smaller values usually lead to a higher sensitivity in that more edges are included into the network. \cr Note that for EBICglasso, an additional specification of the \code{thresh} argument is obsolete, as it is not used for the application of the method.}
#' }
#' (ii) To quantify differences between two networks, the following (a) overall, (b) edge-specific and (c) node-specific network difference characteristics, which have to be supplied in the \code{score.funct} argument, are currently supported:
#' \itemize{
#' \item{\code{frobenius.metric} (overall)} \cr {Calculates the Frobenius metric between two networks}
#' \item{\code{global.str} (overall)} \cr {Calculates the difference in global strength between two networks}
#' \item{\code{maximum.metric} (overall)} \cr {Calculates the maximum metric between two networks}
#' \item{\code{number.differences} (overall)} \cr {Calculates the differences in numbers of edges, clusters and isolated nodes between two networks}
#' \item{\code{spec.dist} (overall)} \cr {Calculates the spectral distance between two networks}
#' \item{\code{jaccard.dist} (overall)} \cr {Calculates the Jaccard distance between two networks}
#' \item{\code{betweenness.inv} (node-specific)} \cr {Calclulates the differences in betweenness between two networks for each node}
#' \item{\code{closeness.inv} (node-specific)} \cr {Calculates the differences in closeness between two networks for each node}
#' \item{\code{degree.inv} (node-specific)} \cr {Calculates the differences in degree between two networks for each node}
#' \item{\code{eigen.inv} (node-specific)} \cr {Calculates the differences in eigenvector centrality between two networks for each node}
#' \item{\code{edge.inv} (edge-specific)} \cr {Calculates the differences in absolute edge weights between two networks for each edge}
#' \item{\code{edge.inv.direc} (edge-specific)} \cr {Calculates the differences in edge weights between two networks for each edge}
#' }
#' Typically, a large number of permutations (e.g. 1000 or 10000) should be chosen in order to obtain reliable results. Note that a large number of permutations may lead to increasing computation times.
#'
#' Note that when underlying data is paired (\code{paired=TRUE}), the input data tables need to have exactly the same dimension, i.e. the same number of columns (nodes) AND rows (samples).
#' @return a list, whose specific structure depends on the specified network difference characteristic in the argument \code{score.funct} (with \eqn{N} denoting the number of nodes):
#' \itemize{
#' \item{if \code{score.funct} is one of \code{frobenius.metric}, \code{global.str}, \code{maximum.metric}, \code{spec.dist} or \code{jaccard.dist}:} \cr {a list with 6 elements: the adjacency matrix for input dat set \code{A} (\eqn{N \times N} matrix), the adjacency matrix for input dat set \code{B} (\eqn{N \times N} matrix), the value of the test statistic (vector of length 1), the values of the test statistics when applying the permutations (vector of length \code{permnum}), the p-value (vector of length 1), the p-value when inserting a pseudocount in the p-value calculation to avoid p-values that are exactly zero in permutation-based settings (vector of length 1)}
#' \item{if \code{score.funct} is \code{number.differences}:} \cr {a list with 6 elements: output when applying \code{create.graph} to input data table \code{A} (list with 15 elements; see documentation of \code{create.graph} function for details), output when applying \code{create.graph} to input data table \code{B} (list with 15 elements; see documentation of \code{create.graph} function for details), the value of the test statistics (vector of length 6), the values of the test statistics when applying the permutations (\code{permnum}\eqn{\times 6} matrix), the p-values (vector of length 6), the p-values when inserting a pseudocount in the p-value calculation to avoid p-values that are exactly zero in permutation-based settings (vector of length 6)}
#' \item{if \code{score.funct} is one of \code{betweenness.inv}, \code{closeness.inv}, \code{degree.inv} or \code{eigen.inv}:} \cr {a list with 6 elements: the adjacency matrix for input dat set \code{A} (\eqn{N \times N} matrix), the adjacency matrix for input dat set \code{B} (\eqn{N \times N} matrix), the value of the test statistic for each node (vector of length \eqn{N}), the values of the test statistics for each node when applying the permutations (\code{permnum}\eqn{\times N} matrix), the p-value for each node (vector of length \eqn{N}), the p-value for each node when inserting a pseudocount in the p-value calculation to avoid p-values that are exactly zero in permutation-based settings (vector of length \eqn{N})}
#' \item{if \code{score.funct} is one of \code{edge.inv} or \code{edge.inv.direc}:} \cr {a list with 8 elements:  the adjacency matrix for input dat set \code{A} (\eqn{N \times N} matrix), the adjacency matrix for input dat set \code{B} (\eqn{N \times N} matrix), the value of the test statistic for each node-node pair (\eqn{N \times N} matrix), the values of the test statistics for each node-node pair when applying the permutations (\code{permnum}\eqn{\times N \times N} array), the p-value for each node-node pair (\eqn{N \times N} matrix), the p-value for each node-node pair when inserting a pseudocount in the p-value calculation to avoid p-values that are exactly zero in permutation-based settings (\eqn{N \times N} matrix), a simplified overview of the p-values for the node-node pairs (\eqn{\frac{N(N-1)}{2} \times 3} matrix; node-node pair in columns 1 and 2, corresponding p-value in column 3), a simplified overview of the p-values for the node-node pairs when inserting a pseudocount in the p-value calculation (\eqn{\frac{N(N-1)}{2} \times 3} matrix; node-node pair in columns 1 and 2, corresponding p-value in column 3)}
#' }
#' @examples
#' ##examples using (a) overall network difference characteristics
#' res1<-perm.test.nw(A=ExDataA,B=ExDataB,permnum=10000,methodlist=list("PCSpearman"),
#' score.funct=frobenius.metric)
#' res2<-perm.test.nw(A=ExDataA,B=ExDataB,permnum=10000,methodlist=list("Spearman"),
#' score.funct=global.str,paired=TRUE)
#'
#' ##examples using (b) node-specific network difference characteristics
#' res3<-perm.test.nw(A=ExDataA,B=ExDataB,permnum=10000,methodlist=list("Spearman"),
#' thresh=0.1,score.funct=betweenness.inv,paired=TRUE)
#' res4<-perm.test.nw(A=ExDataA,B=ExDataB,permnum=10000,methodlist=list("EBICglasso",
#' "spearman",0.1), score.funct=degree.inv)
#'
#' ##examples using (c) edge-specific network difference characteristics
#' res5<-perm.test.nw(A=ExDataA,B=ExDataB,permnum=10000,methodlist=list("Spearman.adj",
#' "bonferroni"),score.funct=edge.inv)
#' res6<-perm.test.nw(A=ExDataA,B=ExDataB,permnum=10000,methodlist=list("EBICglasso",
#' "spearman",0.01),score.funct=edge.inv.direc,paired=TRUE)
#'
#' @export
#'
perm.test.nw<-function (A, B, permnum, methodlist, thresh=0.05, score.funct, paired = FALSE)
{
  stopifnot(`A and B need to be data frames, arrays, matrices or tibbles with the same number of columns (nodes) and column names if applicable.` = all((class(A) %in%  c("matrix", "array","tbl_df",  "tbl", "data.frame") & class(B) %in% c("matrix","array","tbl_df", "tbl", "data.frame") & ncol(A) == ncol(B))), `permnum needs to be a natural number.` = all(permnum >
              0 & permnum%%1 == 0), `thresh needs to be a non-negative number.` = thresh >=
              0, `score.funct needs to be a function.` = class(score.funct) ==
              "function", `paired needs to be a boolean.` = any(paired ==
              TRUE | paired == FALSE))
  if (identical(score.funct, number.differences)) {
    graph.A <- create.graph(A, methodlist = methodlist, thresh = thresh)
    graph.B <- create.graph(B, methodlist = methodlist, thresh = thresh)
    value.orig <- number.differences(graph.A, graph.B)
    names(value.orig) <- c("Graph Diff num of edges",
                           "Graph Diff num of clusters", "Graph Diff num of isolated nodes",
                           "MST Diff num of edges", "MST Diff num of clusters",
                           "MST Diff num of isolated nodes")
  }
  else {
    adj.A <- create.adjacency.matrix(A, methodlist = methodlist,
                                     thresh = thresh)
    adj.B <- create.adjacency.matrix(B, methodlist = methodlist,
                                     thresh = thresh)
    value.orig <- score.funct(adj.A, adj.B)
  }
  if (paired == FALSE) {
    D <- rbind(A, B)
    create.permut.unpaired <- function(seedex) {
      set.seed(seedex)
      s <- sample(nrow(D), replace = FALSE)
      D.new <- D[s, ]
      A.new <- D.new[1:nrow(A), ]
      B.new <- D.new[(nrow(A) + 1):(nrow(D.new)), ]
      output <- list(A.new, B.new)
      return(output)
    }
    shuffle <- lapply(1:permnum, create.permut.unpaired)
  }
  if (paired == TRUE) {
    create.permut.paired <- function(seedex) {
      set.seed(seedex)
      s <- sample(c("A", "B"), nrow(A), replace = TRUE)
      A.new <- rbind(A[s == "A", ], B[s == "B",
      ])
      B.new <- rbind(B[s == "A", ], A[s == "B",
      ])
      output <- list(A.new, B.new)
      return(output)
    }
    shuffle <- lapply(1:permnum, create.permut.paired)
  }
  if (any(identical(score.funct, frobenius.metric), identical(score.funct,
                                                              maximum.metric), identical(score.funct, spec.dist), identical(score.funct,
                                                                                                                            global.str),identical(score.funct, jaccard.dist))) {
    value.perm <- rep(NA, permnum)
    for (n in 1:permnum) {
      adj.A.perm <- create.adjacency.matrix(shuffle[[n]][[1]],
                                            methodlist = methodlist, thresh = thresh)
      adj.B.perm <- create.adjacency.matrix(shuffle[[n]][[2]],
                                            methodlist = methodlist, thresh = thresh)
      value.perm[n] <- score.funct(adj.A.perm, adj.B.perm)
    }
    num.extr <- sum(value.perm >= value.orig)
    pvalue.ecdf <- num.extr/permnum
    pvalue.ecdf.pseudo <- (1 + num.extr)/(permnum + 1)
    output <- list(adj.A, adj.B, value.orig, value.perm,
                   pvalue.ecdf, pvalue.ecdf.pseudo)
    names(output) <- c("adjancency.matrix.A", "adjacency.matrix.B",
                       "test.statistic", "test.statistics.perm",
                       "pvalue", "pvalue.pseudocount")
  }
  if (identical(score.funct, number.differences)) {
    value.perm <- array(data = NA, dim = c(permnum, length(value.orig)))
    for (n in 1:permnum) {
      graph.A.perm <- create.graph(shuffle[[n]][[1]], methodlist = methodlist,
                                   thresh = thresh)
      graph.B.perm <- create.graph(shuffle[[n]][[2]], methodlist = methodlist,
                                   thresh = thresh)
      value.perm[n, ] <- number.differences(graph.A.perm,
                                            graph.B.perm)
    }
    pvalue.ecdf <- rep(NA, length(value.orig))
    pvalue.ecdf.pseudo <- rep(NA, length(value.orig))
    for (j in 1:length(value.orig)) {
      if (all(is.na(value.orig[j]) == TRUE)) {
        pvalue.ecdf[j] <- NA
        pvalue.ecdf.pseudo[j] <- NA
      }
      else {
        if (all(is.na(value.perm[, j]) == FALSE)) {
          num.extr <- sum(value.perm[, j] >= value.orig[j])
          pvalue.ecdf[j] <- num.extr/permnum
          pvalue.ecdf.pseudo[j] <- (1 + num.extr)/(permnum +
                                                     1)
        }
        else {
          permnum <- length(na.omit(value.perm[, j]))
          num.extr <- sum(na.omit(value.perm[, j]) >=
                            value.orig[j])
          pvalue.ecdf[j] <- num.extr/permnum
          pvalue.ecdf.pseudo[j] <- (1 + num.extr)/(permnum +
                                                     1)
        }
      }
    }
    names(pvalue.ecdf) <- c("Graph Diff num of edges",
                            "Graph Diff num of clusters", "Graph Diff num of isolated nodes",
                            "MST Diff num of edges", "MST Diff num of clusters",
                            "MST Diff num of isolated nodes")
    names(pvalue.ecdf.pseudo) <- c("Graph Diff num of edges",
                                   "Graph Diff num of clusters", "Graph Diff num of isolated nodes",
                                   "MST Diff num of edges", "MST Diff num of clusters",
                                   "MST Diff num of isolated nodes")
    output <- list(graph.A, graph.B, value.orig, value.perm,
                   pvalue.ecdf, pvalue.ecdf.pseudo)
    names(output) <- c("graph.A", "graph.B",
                       "test.statistic", "test.statistics.perm",
                       "pvalue", "pvalue.pseudocount")
  }
  if (any(identical(score.funct, degree.inv), identical(score.funct,
                                                        betweenness.inv), identical(score.funct, closeness.inv),
          identical(score.funct, eigen.inv))) {
    value.perm <- array(data = NA, dim = c(permnum, length(value.orig)))
    for (n in 1:permnum) {
      adj.A.perm <- create.adjacency.matrix(shuffle[[n]][[1]],
                                            methodlist = methodlist, thresh = thresh)
      adj.B.perm <- create.adjacency.matrix(shuffle[[n]][[2]],
                                            methodlist = methodlist, thresh = thresh)
      value.perm[n, ] <- score.funct(adj.A.perm, adj.B.perm)
    }
    pvalue.ecdf <- rep(NA, length(value.orig))
    pvalue.ecdf.pseudo <- rep(NA, length(value.orig))
    names(pvalue.ecdf) <- names(value.orig)
    names(pvalue.ecdf.pseudo) <- names(value.orig)
    for (i in 1:length(value.orig)) {
      num.extr <- sum(value.perm[, i] >= value.orig[i])
      pvalue.ecdf[i] <- num.extr/permnum
      pvalue.ecdf.pseudo[i] <- (1 + num.extr)/(permnum +
                                                 1)
    }
    output <- list(adj.A, adj.B, value.orig, value.perm,
                   pvalue.ecdf, pvalue.ecdf.pseudo)
    names(output) <- c("adjacency.matrix.A", "adjacency.matrix.B",
                       "test.statistic", "test.statistics.perm",
                       "pvalue", "pvalue.pseudocount")
  }
  if (any(identical(score.funct, edge.inv), identical(score.funct,
                                                      edge.inv.direc))) {
    value.perm <- array(data = NA, dim = c(permnum, dim(adj.A)))
    for (n in 1:permnum) {
      adj.A.perm <- create.adjacency.matrix(shuffle[[n]][[1]],
                                            methodlist = methodlist, thresh = thresh)
      adj.B.perm <- create.adjacency.matrix(shuffle[[n]][[2]],
                                            methodlist = methodlist, thresh = thresh)
      value.perm[n, , ] <- score.funct(adj.A.perm, adj.B.perm)
    }
    pvalue.ecdf <- array(data = NA, dim = dim(value.orig),
                         dimnames = list(rownames(value.orig), colnames(value.orig)))
    pvalue.ecdf.pseudo <- array(data = NA, dim = dim(value.orig),
                                dimnames = list(rownames(value.orig), colnames(value.orig)))
    for (i in 1:dim(value.orig)[1]) {
      for (j in 1:dim(value.orig)[2]) {
        num.extr <- sum(value.perm[, i, j] >= value.orig[i,
                                                         j])
        pvalue.ecdf[i, j] <- num.extr/permnum
        pvalue.ecdf.pseudo[i, j] <- (1 + num.extr)/(permnum +
                                                      1)
      }
    }
    p.ecdf <- data.frame(var1 = rownames(pvalue.ecdf)[row(pvalue.ecdf)[upper.tri(pvalue.ecdf)]],
                         var2 = colnames(pvalue.ecdf)[col(pvalue.ecdf)[upper.tri(pvalue.ecdf)]],
                         pvalue = pvalue.ecdf[upper.tri(pvalue.ecdf)])
    p.ecdf.pseudo <- data.frame(var1 = rownames(pvalue.ecdf.pseudo)[row(pvalue.ecdf.pseudo)[upper.tri(pvalue.ecdf.pseudo)]],
                                var2 = colnames(pvalue.ecdf.pseudo)[col(pvalue.ecdf.pseudo)[upper.tri(pvalue.ecdf.pseudo)]],
                                pvalue = pvalue.ecdf.pseudo[upper.tri(pvalue.ecdf.pseudo)])
    output <- list(adj.A, adj.B, value.orig, value.perm,
                   pvalue.ecdf, pvalue.ecdf.pseudo, p.ecdf, p.ecdf.pseudo)
    names(output) <- c("adjancency.matrix.A", "adjacency.matrix.B",
                       "test.statistic", "test.statistics.perm",
                       "pvalue", "pvalue.pseudocount", "pvalue.summ",
                       "pvalue.pseudocount.summ")
  }
  return(output)
}

