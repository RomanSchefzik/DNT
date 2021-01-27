#' Plots of the graphs
#' @description
#' This method plots two networks to compare their differences.
#' @param A a matrix.
#' @param B a matrix.
#' @param methodlist the methods used to create the adjacency matrices.
#' @param thresh the threshold under which the values of the correlations are relevant.
#' @param networkA Boolean value that specifies if the Survivors Admission shall be plotted.
#' @param networkB Boolean value that specifies if the Survivors Event shall be plotted.
#' @param networkAtitle The title of the graph of the network of A.
#' @param networkBtitle The title of the graph of the network of B.
#' @param cluster is the boolean value that states if the graphs shall be divided in their clusters or not.
#' @param negcol the color of the edges with a negative correlation value.
#' @param poscol the color of the edges with a positive correlation value.
#' @param multiplier the number multiplied with the edge weight to get a bigger differnece between the thick and thin lines.
#' @param layout the layout of the networks.
#' @param vSize the size of the vertices.
#' @param tSize the size of the names inside the vertices.
#' @details
#' The graphs are plotted with the same layout and other plot characteristics so that the two graphs are easy to
#' compare. The left network is the Survivors Admission and the right one the Survivor Event.
#' @return
#' @examples
#' graph.plot(x1,x2, list("Spearman), caption = "Spearman", directory)
#' @export
#'
graph.plot <- function(A, B, methodlist, thresh = 0.05, networkA = TRUE, networkB = TRUE, networkAtitle = "Network A", networkBtitle = "Network B", cluster = TRUE, negcol = "red", poscol= "blue",multiplier = 4, curved = TRUE, layout = layout.fruchterman.reingold, vSize = 16, tSize = 0.8){

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

#' Create cluster of an Igraph
#' @description
#' The method creates an Igraph out of the matrix A and the Igraphs clusters.
#' @param A a matrix.
#' @param methodlist the methods used to create the adjacency matrix.
#' @param thresh the threshold under which the values of the correlations are relevant.
#' @details
#' First the method creates the adjacency matrix out of the matrix A and after that the adjacency matrix can be
#' used to create the Igraph and its clusters.
#' @return The method returns a list with the Igraph from A and its clusters.
#' @example
#' create.Igraphclustering(x1, list("Spearman"))
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
  }
  list(g.clustering.x1, g.x1)
}
