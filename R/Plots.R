#' Plots of the graphs
#' @description
#' This method plots two networks to compare their differences.
#' @param A a matrix.
#' @param B a matrix.
#' @param methodlist the methods used to create the adjacency matrices.
#' @param thresh the threshold under which the values of the correlations are relevant.
#' @param admission Boolean value that specifies if the Survivors Admission shall be plotted.
#' @param event Boolean value that specifies if the Survivors Event shall be plotted.
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
graph.plot <- function(A, B, methodlist, thresh = 0.05, admission = TRUE, event = TRUE, cluster = TRUE, negcol = "red", poscol= "blue",multiplier = 4, curved = TRUE, layout = layout.fruchterman.reingold, vSize = 16, tSize = 0.8){

  stopifnot("cluster needs to be boolean" = cluster==TRUE || cluster==FALSE,
            "negcol and poscol need to be strings" = class(negcol)=="character" && class(poscol)=="character",
            "multiplier needs to be a positive number" = (class(multiplier)=="numeric"||class(multiplier)=="integer") && multiplier>0,
            "admission and event needs to be boolean" = class(admission)=="logical" && class(event)=="logical",
            "layout needs to be a layout-function" = class(layout)=="function",
            "vSize and tSize need to be positive numbers" = (class(vSize)=="numeric" || class(vSize)=="integer") &&
                                                            (class(tSize)=="numeric" || class(tSize)=="integer"),
            "vSize needs to be 20 times bigger than tSize to have a good relation between the names and the vertices" = tSize*20<=vSize
            )

  g1 <- create.Igraphclustering(A, methodlist, thresh)
  g2 <- create.Igraphclustering(B, methodlist, thresh)

  if(is.null(g1[[1]]) == FALSE){
    E(g1[[2]])[which(E(g1[[2]])$weight<0)]$color <- negcol
    E(g1[[2]])[which(E(g1[[2]])$weight>0)]$color <- poscol

    E(g1[[2]])$weight <- abs(E(g1[[2]])$weight)
  }

  if(is.null(g2[[1]]) == FALSE){
    E(g2[[2]])[which(E(g2[[2]])$weight<0)]$color <- negcol
    E(g2[[2]])[which(E(g2[[2]])$weight>0)]$color <- poscol

    E(g2[[2]])$weight <- abs(E(g2[[2]])$weight)
  }

  par(mfrow=c(1,2),oma=c(0,0,2,0))


  if(cluster == TRUE){
    if(admission == TRUE && is.null(g1[[1]]) == FALSE){
      plot(
        g1[[1]], g1[[2]],
        layout= layout,
        edge.curved=curved,
        vertex.size=vSize,
        vertex.label.cex=tSize,
        edge.width=E(g1[[2]])$weight*multiplier,
        edge.color = E(g1[[2]])$color,
        edge.arrow.mode=FALSE,
        main="Survivors Admission")
    }else if(admission == TRUE){
      plot(
        g1[[2]],
        layout=layout,
        edge.curved=curved,
        vertex.size=vSize,
        vertex.label.cex=tSize,
        edge.width=E(g1[[2]])$weight*multiplier,
        edge.arrow.mode=FALSE,
        main="Survivors Admission")
      mtext("No correlations!", side = 1)
    }

    if(event == TRUE && is.null(g2[[1]]) == FALSE){
      plot(
        g2[[1]], g2[[2]],
        layout=layout,
        edge.curved=curved,
        vertex.size=vSize,
        vertex.label.cex=tSize,
        edge.width=E(g2[[2]])$weight*multiplier,
        edge.color = E(g2[[2]])$color,
        edge.arrow.mode=FALSE,
        main="Survivors Event")
    }else if(event == TRUE){
      plot(
        g2[[2]],
        layout=layout,
        edge.curved=curved,
        vertex.size=vSize,
        vertex.label.cex=tSize,
        edge.width=E(g2[[2]])$weight*multiplier,
        edge.arrow.mode=FALSE,
        main="Survivors Event")
      mtext("No correlations!", side = 1)
    }

  }else{
    if(admission == TRUE){
      plot(
        g1[[2]],
        layout=layout,
        edge.curved=curved,
        vertex.size=vSize,
        vertex.label.cex=tSize,
        edge.width=E(g1[[2]])$weight*multiplier,
        edge.arrow.mode=FALSE,
        main="Survivors Admission")
    }

    if(event == TRUE){
      plot(
        g2[[2]],
        layout=layout,
        edge.curved=curved,
        vertex.size=vSize,
        vertex.label.cex=tSize,
        edge.width=E(g2[[2]])$weight*multiplier,
        edge.arrow.mode=FALSE,
        main="Survivors Event")
    }
  }
  caption <- methodlist[[1]]
  mtext(caption,outer=TRUE,cex=1.5)

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
create.Igraphclustering <- function(A, methodlist, thresh = 0.05){

  cm.x1 <- create.adjacency.matrix(A, methodlist, thresh)

  if(all(cm.x1==0) == TRUE){
    numbervertices <- length(A)
    g.x1 <- make_empty_graph(n = numbervertices)
    V(g.x1)$name <- names(A)
    g.clustering.x1 <- NULL
  }
  else{
    # Make an Igraph object from this matrix:
    g.x1<-graph.adjacency(cm.x1,weighted=TRUE, mode="undirected", diag=FALSE)

    # Simplfy the adjacency object
    g.x1<-simplify(g.x1, remove.multiple=TRUE, remove.loops=TRUE)

    # let's see if we have communities here using the
    # Grivan-Newman algorithm
    # 1st we calculate the edge betweenness, merges, etc...
    g.communities.x1 <- edge.betweenness.community(g.x1, weights=NULL, directed=FALSE)
    g.clustering.x1 <- make_clusters(g.x1, membership=g.communities.x1$membership)
    V(g.x1)$color <- g.communities.x1$membership
  }
  list(g.clustering.x1, g.x1)
}
