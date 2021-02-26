#' Provides plots allowing for visual comparison of two networks
#' @description
#' This function provides plots allowing for visual comparison of two networks, which are given by adjacency matrices that are produced out of input data tables, using an estimation method specified by the user.
#' @param A,B input data tables from which the adjacency matrices will be generated, to be provided in form of matrices, arrays, data frames or tibbles
#' @param methodlist a list specifying the method which is used to estimate and create the adjacency matrix; see details for further information
#' @param thresh a number between 0 and 1 (default is set to 0.05) specifying the singificance level: if the p-value corresponding to an edge weight is greater than \code{thresh}, the corresponding edge weight is not considered to be significant and thus set to zero
#' @param networkA Boolean, indicating whether the network corresponding to the first data set shall be plotted
#' @param networkB Boolean, indicating whether the network corresponding to the second data set shall be plotted
#' @param networkAtitle a character, indicating the title of the first network
#' @param networkBtitle a character, indicating the title of the second network
#' @param cluster a Boolean, indicating whether the networks shall be plotted along with the derived communities/clusters
#' @param negcol a character, specifying the color which is used to represent edges with a negative correlation
#' @param poscol a character, specifying the color which is used to represent edges with a positive correlation
#' @param multiplier a number, representing the factor with which an edge weight is multiplied in order to regulate the thickness of the drawn edges
#' @param curved a Boolean, indicating whether the edges should be drawn using a curved or solid line
#' @param layout a function specifying the layout of the networks; see details for possible options and further information
#' @param vSize a number, specifying the size of the nodes
#' @param tSize a number, specifying the node label size
#' @details
#' This function provides plots allowing for visual comparison of two networks, which are given by adjacency matrices that are produced out of input data tables, using an estimation method specified by the user. The network estimation method has to be specified in form of a list in the \code{methodlist} argument. Currently, the following estimation methods are supported:
#' \itemize{
#' \item{\code{list("Spearman")} \cr Edge weights are estimated using Spearman correlation, where unadjusted p-values are employed to determine significance. To apply this method, the expression \code{list("Spearman")} has to be provided in the \code{methodlist} argument.}
#' \item{\code{list("Spearman.adj",adjustment method)} \cr Edge weights are estimated using Spearman correlation, where p-values adjusted for multiple testing are employed to determine significance. To apply this method, the expression \code{list("Spearman.adj",adjustment method)} has to be provided in the \code{methodlist} argument, where \code{adjustment method} has to be one of the options for multple testing adjustment provided by the standard \code{p.adjust} R function, i.e. one of \code{"BH"}, \code{"bonferroni"}, \code{"BY"}, \code{"fdr"}, \code{"hochberg"}, \code{"holm"} or \code{"hommel"}.}
#' \item{\code{list("PCSpearman")} \cr Edge weights are estimated using partial Spearman correlation, where unadjusted p-values are employed to determine significance. To apply this method, the expression \code{list("PCSpearman")} has to be provided in the \code{methodlist} argument.}
#' \item{\code{list("PCSpearman.adj",adjustment method)} \cr Edge weights are estimated using partial Spearman correlation, where p-values adjusted for multiple testing are employed to determine significance. To apply this method, the expression \code{list("PCSpearman.adj",adjustment method)} has to be provided in the \code{methodlist} argument, where \code{adjustment method} has to be one of the options for multple testing adjustment provided by the standard \code{p.adjust} R function, i.e. one of \code{"BH"}, \code{"bonferroni"}, \code{"BY"}, \code{"fdr"}, \code{"hochberg"}, \code{"holm"} or \code{"hommel"}.}
#' \item{\code{list("DistCorr")} \cr Edge weights are estimated using distance correlation, where unadjusted p-values are employed to determine significance. To apply this method, the expression \code{list("DistCorr")} has to be provided in the \code{methodlist} argument. Note that the calculations may require larger computation times, as a permutation test is involved to derive the corresponding p-values for the distance correlations.}
#' \item{\code{list("DistCorr.adj",adjustment method)} \cr Edge weights are estimated using distance correlation, where p-values adjusted for multiple testing are employed to determine significance. To apply this method, the expression \code{list("DistCorr.adj",adjustment method)} has to be provided in the \code{methodlist} argument, where \code{adjustment method} has to be one of the options for multple testing adjustment provided by the standard \code{p.adjust} R function, i.e. one of \code{"BH"}, \code{"bonferroni"}, \code{"BY"}, \code{"fdr"}, \code{"hochberg"}, \code{"holm"} or \code{"hommel"}. Note that the calculations may require larger computation times, as a permutation test is involved to derive the corresponding p-values for the distance correlations.}
#' \item{\code{list("EBICglasso",correlation type,tuning parameter)} \cr Edge weights are estimated using the EBICglasso approach. To apply this method, the expression \code{list("EBICglasso",correlation type,tuning parameter)} has to be provided in the \code{methodlist} argument. Here, \code{correlation type} has to be one of the association options provided by the standard \code{cor} R function, i.e. one of \code{"kendall"}, \code{"pearson"} or \code{"spearman"}. Moreover, \code{tuning parameter} has to be a number specifying the EBIC tuning parameter \eqn{\gamma}. Typical choices include values between 0 and 0.5, where smaller values usually lead to a higher sensitivity in that more edges are included into the network. \cr Note that for EBICglasso, an additional specification of the \code{thresh} argument is obsolete, as it is not used for the application of the method.}
#' }
#' The following options (functions) from the \code{igraph} R package are provided to specify the layout of the plots in the \code{layout} argument:
#' \itemize{
#' \item{\code{igraph::layout.auto} (default)}
#' \item{\code{igraph::layout.circle}}
#' \item{\code{igraph::layout.davidson.harel}}
#' \item{\code{igraph::layout.drl}}
#' \item{\code{igraph::layout.fruchterman.reingold}}
#' \item{\code{igraph::layout.gem}}
#' \item{\code{igraph::layout.graphopt}}
#' \item{\code{igraph::layout.grid}}
#' \item{\code{igraph::layout.kamada.kawai}}
#' \item{\code{igraph::layout.lgl}}
#' \item{\code{igraph::layout.mds}}
#' \item{\code{igraph::layout.reingold.tilfort}}
#' \item{\code{igraph::layout.star}}
#' \item{\code{igraph::layout.svd}}
#' }
#' @return a plot of the two specified networks, which allows for visual comparison
#' @examples
#' comp.plot(ExDataA,ExDataB,methodlist=list("Spearman"))
#' comp.plot(ExDataA,ExDataB,methodlist=list("PCSpearman.adj","bonferroni"),
#' layout=igraph::layout.circle,curved=FALSE)
#' comp.plot(ExDataA,ExDataB,methodlist=list("EBICglasso","spearman",0.1),
#' layout=igraph::layout.fruchterman.reingold,curved=FALSE)
#' comp.plot(ExDataA,ExDataB,methodlist=list("EBICglasso","pearson",0.05),
#' layout=igraph::layout.star,cluster=FALSE)
#' @export
#'
comp.plot<-function(A, B, methodlist, thresh = 0.05, networkA = TRUE, networkB = TRUE, networkAtitle = "Network A", networkBtitle = "Network B", cluster = TRUE, negcol = "red", poscol= "blue",multiplier = 4, curved = TRUE, layout = igraph::layout.auto, vSize = 16, tSize = 0.8){

  stopifnot(`cluster needs to be boolean` = cluster==TRUE || cluster==FALSE,
            `negcol and poscol need to be strings` = class(negcol)=="character" && class(poscol)=="character",
            `multiplier needs to be a positive number` = class(multiplier) %in% c("numeric","integer") && multiplier>0,
            `networkA and networkB must be boolean` = class(networkA)=="logical" && class(networkB)=="logical",
            `networkAtitle and networkBtitle must be srings` = class(networkAtitle)=="character" && class(networkBtitle)=="character",
            `layout needs to be a layout-function` = class(layout)=="function",
            `vSize and tSize need to be positive numbers` = class(vSize) %in% c("numeric","integer") &&
              class(tSize) %in% c("numeric","integer"),
            `vSize needs to be 20 times bigger than tSize to have a good relation between the names and the vertices` = tSize*20<=vSize)

  if(networkA == TRUE){
    g1 <- create.Igraphclustering(A, methodlist, thresh)

    if(is.null(g1[[1]]) == FALSE){
      igraph::E(g1[[2]])[which(igraph::E(g1[[2]])$weight<0)]$color <- negcol
      igraph::E(g1[[2]])[which(igraph::E(g1[[2]])$weight>0)]$color <- poscol

      igraph::E(g1[[2]])$weight <- abs(igraph::E(g1[[2]])$weight)
    }
  }
  if(networkB == TRUE){
    g2 <- create.Igraphclustering(B, methodlist, thresh)

    if(is.null(g2[[1]]) == FALSE){
      igraph::E(g2[[2]])[which(igraph::E(g2[[2]])$weight<0)]$color <- negcol
      igraph::E(g2[[2]])[which(igraph::E(g2[[2]])$weight>0)]$color <- poscol

      igraph::E(g2[[2]])$weight <- abs(igraph::E(g2[[2]])$weight)
    }
  }


  graphics::par(mfrow=c(1,2),oma=c(0,0,2,0))

  if(cluster == TRUE){
    if(networkA == TRUE && is.null(g1[[1]]) == FALSE){
      plot(
        g1[[1]], g1[[2]],
        layout= layout,
        edge.curved=curved,
        vertex.size=vSize,
        vertex.label.cex=tSize,
        edge.width=E(g1[[2]])$weight*multiplier,
        edge.color = E(g1[[2]])$color,
        edge.arrow.mode=FALSE,
        main=networkAtitle)
    }else if(networkA == TRUE){
      plot(
        g1[[2]],
        layout=layout,
        edge.curved=curved,
        vertex.size=vSize,
        vertex.label.cex=tSize,
        edge.width=E(g1[[2]])$weight*multiplier,
        edge.arrow.mode=FALSE,
        main=networkAtitle)
      mtext("No correlations!", side = 1)
    }

    if(networkB == TRUE && is.null(g2[[1]]) == FALSE){
      plot(
        g2[[1]], g2[[2]],
        layout=layout,
        edge.curved=curved,
        vertex.size=vSize,
        vertex.label.cex=tSize,
        edge.width=E(g2[[2]])$weight*multiplier,
        edge.color = E(g2[[2]])$color,
        edge.arrow.mode=FALSE,
        main=networkBtitle)
    }else if(networkB == TRUE){
      plot(
        g2[[2]],
        layout=layout,
        edge.curved=curved,
        vertex.size=vSize,
        vertex.label.cex=tSize,
        edge.width=E(g2[[2]])$weight*multiplier,
        edge.arrow.mode=FALSE,
        main=networkBtitle)
      mtext("No correlations!", side = 1)
    }

  }else{
    if(networkA == TRUE){
      plot(
        g1[[2]],
        layout=layout,
        edge.curved=curved,
        vertex.size=vSize,
        vertex.label.cex=tSize,
        edge.width=E(g1[[2]])$weight*multiplier,
        edge.arrow.mode=FALSE,
        main=networkAtitle)
    }

    if(networkB == TRUE){
      plot(
        g2[[2]],
        layout=layout,
        edge.curved=curved,
        vertex.size=vSize,
        vertex.label.cex=tSize,
        edge.width=E(g2[[2]])$weight*multiplier,
        edge.arrow.mode=FALSE,
        main=networkBtitle)
    }
  }
  caption <- methodlist[[1]]
  graphics::mtext(caption,outer=TRUE,cex=1.5)
}
utils::globalVariables("plot")

#' Creates \code{igraph} graph and clustering results
#' @description
#' This function creates \code{igraph} graph and clustering results from an adjacency matrix that is produced out of an input data table, using an estimation method specified by the user.
#' @param A input data table from which the adjacency matrix will be generated, to be provided in form of a matrix, array, data frame or tibble
#' @param methodlist a list specifying the method which is used to estimate and create the adjacency matrix; see details for further information
#' @param thresh a number between 0 and 1 (default is set to 0.05) specifying the singificance level: if the p-value corresponding to an edge weight is greater than \code{thresh}, the corresponding edge weight is not considered to be significant and thus set to zero
#' @details
#' This function creates \code{igraph} graph and clustering results from an adjacency matrix that is produced out of an input data table, using an estimation method specified by the user. Clustering is performed using the Girvan-Newman algorithm based on edge betweenness. The function is used in the central \code{comp.plot} function of the DNT package for visual comparison of two networks and builds on respective implementations in the \code{igraph} R package. \cr
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
#' @return a list with two elements: the first list element contains the results for the clustering by the Girvan-Newman algorithm based on edge betweenness; the second list element contains the results for the creation of an \code{igraph} graph
#' @examples
#' create.Igraphclustering(ExDataA,methodlist=list("Spearman"))
#' create.Igraphclustering(ExDataA,methodlist=list("Spearman.adj","bonferroni"))
#' create.Igraphclustering(ExDataA,methodlist=list("PCSpearman"),thresh=0.1)
#' create.Igraphclustering(ExDataA,methodlist=list("PCSpearman.adj","BH"),thresh=0.1)
#' create.Igraphclustering(ExDataA,methodlist=list("DistCorr"))
#' create.Igraphclustering(ExDataA,methodlist=list("DistCorr.adj","bonferroni"))
#' create.Igraphclustering(ExDataB,methodlist=list("EBICglasso","spearman",0.1))
#' create.Igraphclustering(ExDataB,methodlist=list("EBICglasso","pearson",0.05))
#' @export
#'
create.Igraphclustering<-function (A, methodlist, thresh = 0.05)
{
  cm.x1 <- create.adjacency.matrix(A, methodlist, thresh)
  if (all(cm.x1 == 0) == TRUE) {
    numbervertices <- length(A)
    g.x1 <- igraph::make_empty_graph(n = numbervertices)
    igraph::V(g.x1)$name <- names(A)
    g.clustering.x1 <- NULL
  }
  else {
    g.x1 <- igraph::graph.adjacency(cm.x1, weighted = TRUE, mode = "undirected",
                            diag = FALSE)
    g.x1 <- igraph::simplify(g.x1, remove.multiple = TRUE, remove.loops = TRUE)
    g.communities.x1 <- igraph::edge.betweenness.community(g.x1,
                                                   weights = NULL, directed = FALSE)
    g.clustering.x1 <- igraph::make_clusters(g.x1, membership = g.communities.x1$membership)
    igraph::V(g.x1)$color <- g.communities.x1$membership
    igraph::V(g.x1)$label.color <- "black"
  }
  list(g.clustering.x1, g.x1)
}

