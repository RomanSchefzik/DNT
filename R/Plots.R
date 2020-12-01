#' Plots of the graphs
#' @description
#' This method plots two networks to compare their differences.
#' @param A a matrix.
#' @param B a matrix.
#' @param methodlist the methods used to create the adjacency matrices.
#' @param cluster is the boolean value that states if the graphs shall be divided in their clusters or not.
#' @param negcol the color of the edges with a negative correlation value.
#' @param poscol the color of the edges with a positive correlation value.
#' @param caption the title of the pdf-file.
#' @param layout the layout of the networks.
#' @param vSize the size of the vertices.
#' @param tSize the size of the names inside the vertices.
#' @param directory the directory where the file shall be saved.
#' @details
#' The graphs are plotted with the same layout and other plot characteristics so that the two graphs are easy to
#' compare. The left network is the Survivors Admission and the right one the Survivor Event.
#' @return
#' @examples
#' graph.plot(x1,x2, list("Spearman), caption = "Spearman", directory)
#' @export
#'
graph.plot <- function(A, B, methodlist, cluster = TRUE, negcol = "red", poscol= "blue", caption, layout = layout.fruchterman.reingold, vSize = 16, tSize = 0.8, directory){

  stopifnot("cluster needs be boolean" = cluster==TRUE || cluster==FALSE,
            "negcol and poscol need to be strings" = class(negcol)=="character" && class(poscol)=="character",
            "caption needs to be a string" = class(caption)=="character",
            "layout needs to be a layout-function" = class(layout)=="function",
            "vSize and tSize need to be positive numbers" = (class(vSize)=="numeric" || class(vSize)=="integer") &&
                                                            (class(tSize)=="numeric" || class(tSize)=="integer"),
            "vSize needs to be 20 times bigger than tSize to have a good relation between the names and the vertices" = tSize*20<=vSize,
            "directory needs to be a string" = class(directory)=="character")

  g1 <- create.Igraphclustering(A, methodlist)
  g2 <- create.Igraphclustering(B, methodlist)

  E(g1[[2]])[which(E(g1[[2]])$weight<0)]$color <- negcol
  E(g1[[2]])[which(E(g1[[2]])$weight>0)]$color <- poscol

  E(g2[[2]])[which(E(g2[[2]])$weight<0)]$color <- negcol
  E(g2[[2]])[which(E(g2[[2]])$weight>0)]$color <- poscol

  E(g1[[2]])$weight <- abs(E(g1[[2]])$weight)
  E(g2[[2]])$weight <- abs(E(g2[[2]])$weight)

  filename <- stringr::str_glue(directory, ".pdf", sep ="")
  pdf(file=filename,width=15,height=8)#,pointsize=15)

  par(mfrow=c(1,2),oma=c(0,0,2,0))


  if(cluster == TRUE){
    plot(
      g1[[1]], g1[[2]],
      layout= layout,
      edge.curved=TRUE,
      vertex.size=vSize,
      #vertex.label.dist=-0.5,
      #vertex.label.color="red",
      #asp=1,
      vertex.label.cex=tSize,
      edge.width=E(g1[[2]])$weight*4,
      edge.color = E(g1[[2]])$color,
      edge.arrow.mode=FALSE,
      main="Survivors Admission")


    plot(
      g2[[1]], g2[[2]],
      layout=layout,
      edge.curved=TRUE,
      vertex.size=vSize,
      #vertex.label.dist=-0.5,
      #vertex.label.color="black",
      #asp=FALSE,
      vertex.label.cex=tSize,
      edge.width=E(g2[[2]])$weight*4,
      edge.color = E(g2[[2]])$color,
      edge.arrow.mode=FALSE,
      main="Survivors Event")

  }else{
    plot(
      g1[[2]],
      layout=layout,
      edge.curved=TRUE,
      vertex.size=vSize,
      #vertex.label.dist=-0.5,
      #vertex.label.color=E(g1[[2]])$color,
      #asp=FALSE,
      vertex.label.cex=tSize,
      edge.width=E(g1[[2]])$weight*4,
      edge.arrow.mode=FALSE,
      main="Survivors Admission")


    plot(
      g2[[2]],
      layout=layout,
      edge.curved=TRUE,
      vertex.size=vSize,
      #vertex.label.dist=-0.5,
      #vertex.label.color="black",
      #asp=FALSE,
      vertex.label.cex=tSize,
      edge.width=E(g2[[2]])$weight*4,
      edge.arrow.mode=FALSE,
      main="Survivors Event")

  }

  #caption <- methodlist[[1]]
  mtext(caption,outer=TRUE,cex=1.5)

  dev.off()
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

  # Make an Igraph object from this matrix:
  g.x1<-graph.adjacency(cm.x1,weighted=TRUE, mode="undirected", diag=FALSE)

  # Simplfy the adjacency object
  g.x1<-simplify(g.x1, remove.multiple=TRUE, remove.loops=TRUE)

  # Colour negative correlation edges as red
  #E(g.x1)[which(E(g.x1)$weight<0)]$color <- "red"

  # Colour positive correlation edges as blue
  #E(g.x1)[which(E(g.x1)$weight>0)]$color <- "blue"

  # Convert edge weights to absolute values
  #E(g.x1)$weight <- abs(E(g.x1)$weight)

  # let's see if we have communities here using the
  # Grivan-Newman algorithm
  # 1st we calculate the edge betweenness, merges, etc...
  g.communities.x1 <- edge.betweenness.community(g.x1, weights=NULL, directed=FALSE)
  g.clustering.x1 <- make_clusters(g.x1, membership=g.communities.x1$membership)
  V(g.x1)$color <- g.communities.x1$membership

  list(g.clustering.x1, g.x1)
}
