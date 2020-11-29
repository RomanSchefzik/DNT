
###1.) Network estimation

#' Create Adjacency Matrix
#' @description
#' The method creates an adjacency matrix out of the matrix A with the methods specified in the attribute methodlist.
#' @param A a matrix.
#' @param methodlist the methods used to create the adjacency matrix.
#' @param thresh the threshold under which the values of the correlations are relevant.
#' @details
#' To show the meaning of a certain value in the adjacency matrix, the titles of the columns of the matrix A are used
#' as titles of the rows and the columns of the adjacency matrix.
#' @return The adjacency matrix resulting from the matrix A.
#' @examples
#' create.adjacency.matrix(A,list("Spearman"))
#' create.adjacency.matrix(A,list("DistCorr.adj", "holm"))
#' @export
#'
create.adjacency.matrix<-function(A,methodlist,thresh=0.05){

  stopifnot("A needs to be a tibble or a matrix" = class(A)==c("tbl_df", "tbl", "data.frame")||class(A)==c("matrix","array"),
            "The matrix A needs to have more than 4 rows" = nrow(A)>4,
            "The matrix A needs to have at least one column" = ncol(A)>=1,
            "Methodlist needs to be a list of Strings  with at least one element (except for the method EBICglasso where the third elements needs to be a number)" =
                    class(methodlist)=="list" && length(methodlist)>=1 &&
                    (all(lapply(methodlist,class)=="character") || methodlist[[1]]=="EBICglasso" && class(methodlist[[2]])=="character" && (class(methodlist[[3]])=="numeric"||class(methodlist[[3]])=="integer")),
            "For methods with the ending '.adj' a second method is needed and when using the method 'EBICglasso' one extra method and one extra value are needed" =
                    (stringr::str_detect(methodlist[[1]], "\\.adj") == FALSE && length(methodlist)>=1)
                    || (stringr::str_detect(methodlist[[1]], "\\.adj") && length(methodlist)>=2)
                    || (methodlist[[1]] == "EBICglasso" && length(methodlist)>=3),
            "Thresh needs to be a non-negative number" = class(thresh)=="numeric" && thresh>=0)

  method <- methodlist[[1]]
  if(method=="Spearman"){
    correls<-rcorr(as.matrix(A),type="spearman")
    cm<-correls$r
    # Keep only correlations with p-value <=thresh
    diag(cm)<-0
    cm[correls$P>thresh]<-0
  }

  if(method=="PCSpearman"){
    correls<-pcor(as.matrix(A),method="spearman")
    cm<-correls$estimate
    # Keep only correlations with p-value <=thresh
    diag(cm)<-0
    cm[correls$p.value>thresh]<-0
  }

  if(method=="Spearman.adj"){
    correls<-rcorr(as.matrix(A),type="spearman")
    cm<-correls$r

    pmat<-correls$P
    ltri<-lower.tri(pmat)
    utri<-upper.tri(pmat)

    adj.method <- methodlist[[2]]
    pmat[ltri]<-p.adjust(pmat[ltri], method = adj.method)
    pmat[utri]<-t(pmat)[utri]

    # Keep only correlations with adjusted p-value <=thresh
    diag(cm)<-0
    cm[pmat>thresh]<-0
  }

  if(method=="PCSpearman.adj"){
    correls<-pcor(as.matrix(A),method="spearman")
    cm<-correls$estimate

    pmat<-correls$p.value
    ltri<-lower.tri(pmat)
    utri<-upper.tri(pmat)
    adj.method <- methodlist[[2]]
    pmat[ltri]<-p.adjust(pmat[ltri], method = adj.method)
    pmat[utri]<-t(pmat)[utri]

    # Keep only correlations with adjusted p-value <=thresh
    diag(cm)<-0
    cm[pmat>thresh]<-0
  }


  if(method=="DistCorr"){
    A<-as.matrix(A)

    cm<-array(data=NA,dim=c(dim(A)[2],dim(A)[2]))
    for(i in 1:dim(A)[2]){
      for(j in 1:dim(A)[2]){
        cm[i,j]<-dcor(A[,i],A[,j],index=1.0)
      }}
    rownames(cm)<-colnames(A)
    colnames(cm)<-colnames(A)

    ##compute p-values based on 10000 permutations
    pvals<-array(data=NA,dim=c(dim(A)[2],dim(A)[2]))
    for(i in 1:dim(A)[2]){
      for(j in 1:dim(A)[2]){
        set.seed(24)
        pvals[i,j]<-dcor.test(A[,i],A[,j],index=1.0,R=10000)$p.value
      }}
    rownames(pvals)<-colnames(A)
    colnames(pvals)<-colnames(A)

    # Keep only correlations with p-value <=thresh
    diag(cm)<-0
    cm[pvals>thresh]<-0
  }

  if(method=="DistCorr.adj"){
    A<-as.matrix(A)

    cm<-array(data=NA,dim=c(dim(A)[2],dim(A)[2]))
    for(i in 1:dim(A)[2]){
      for(j in 1:dim(A)[2]){
        cm[i,j]<-dcor(A[,i],A[,j],index=1.0)
      }}
    rownames(cm)<-colnames(A)
    colnames(cm)<-colnames(A)

    ##compute p-values based on 10000 permutations
    pvals<-array(data=NA,dim=c(dim(A)[2],dim(A)[2]))
    for(i in 1:dim(A)[2]){
      for(j in 1:dim(A)[2]){
        set.seed(24)
        pvals[i,j]<-dcor.test(A[,i],A[,j],index=1.0,R=10000)$p.value
      }}
    rownames(pvals)<-colnames(A)
    colnames(pvals)<-colnames(A)

    pmat<-pvals
    ltri<-lower.tri(pmat)
    utri<-upper.tri(pmat)
    adj.method <- methodlist[[2]]
    pmat[ltri]<-p.adjust(pmat[ltri], method = adj.method)
    pmat[utri]<-t(pmat)[utri]
    rownames(pmat)<-colnames(A)
    colnames(pmat)<-colnames(A)

    # Keep only correlations with adjusted p-value <=thresh
    diag(cm)<-0
    cm[pmat>thresh]<-0
  }


  if(method=="EBICglasso"){
    corm.ebicglasso <- methodlist[[2]]
    tun.ebicglasso <- methodlist[[3]]
    correls<-cor(A,method=corm.ebicglasso)
    cm<-EBICglasso(correls,n=dim(A)[1],gamma=tun.ebicglasso)
    # Keep only correlations with p-value <=thresh; NOT NECESSARY FOR EBICGLASSO, AS VARIABLE SELECTION IS PERFORMED
    diag(cm)<-0
    cm[correls$p.value>thresh]<-0
  }

  output<-cm
  return(output)

}


#' Create Graph
#' @description
#' The method creates an graph out of the matrix A.
#' @param A a matrix.
#' @param methodlist the methods used to create the adjacency matrix.
#' @param thresh the threshold under which the values of the correlations are relevant.
#' @details
#' To create the graph (adjacency) out of any matrix A this matrix needs to be converted into an adjacency matrix. This is
#' done by the method create.adjacency.matrix. Then many of the graphs characteristics are evaluated and saved in variables.
#' With the help of Prim's Algorithm the graph (adjacency) is converted into a minimum spanning tree (short: MST) of which
#' the same characteristics are evaluated.
#' @return create.graph returns a list of the adjacency matrix, the graph, the graph communities, the graph clustering and
#' the graphs vertex degrees, number of edges, number of clusters and number of isolated nodes. Followed by the minimum
#' spanning tree, its communities, its clustering and the MSTs vertex degrees, number of edges, number of clusters and
#' number of isolated nodes.
#' @examples
#' create.graph(A,list("Spearman"))
#' @export
#'
create.graph<-function(A,methodlist,thresh=0.05){

  cm <- create.adjacency.matrix(A, methodlist,thresh)

  if(all(cm==0)){
    output<-list(cm,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
    names(output)<-c("adjacency matrix","graph","graph communities","graph clustering","graph vertex degrees","graph number of edges","graph number of clusters","graph number of isolated nodes","MST","MST communities","MST clustering","MST vertex degrees","MST number of edges","MST number of clusters","MST number of isolated nodes")
    return(output)
  }else{

    # Make an Igraph object from this matrix:
    g<-graph.adjacency(cm,weighted=TRUE, mode="undirected", diag=FALSE)

    # Simplfy the adjacency object
    g<-simplify(g, remove.multiple=TRUE, remove.loops=TRUE)

    # Colour negative correlation edges as red
    E(g)[which(E(g)$weight<0)]$color <- "red"

    # Colour positive correlation edges as blue
    E(g)[which(E(g)$weight>0)]$color <- "blue"

    # Convert edge weights to absolute values
    # E(g)$weight <- abs(E(g)$weight)

    # Assign names to the graph vertices (optional)
    V(g)$name <- V(g)$name


    # let's see if we have communities here using the
    # Grivan-Newman algorithm
    # 1st we calculate the edge betweenness, merges, etc...
    g.communities <- edge.betweenness.community(g, weights=NULL, directed=FALSE)
    g.clustering <- make_clusters(g, membership=g.communities$membership)
    V(g)$color <- g.communities$membership

    # Check the vertex degree, i.e., number of connections to each vertex
    deg.g<-degree(g)

    ##number of edges
    edgenum.g<-gsize(g)


    ###number of clusters
    clusnum.g<-length(unique(g.communities$membership))

    ##number of isolated nodes (i.e. nodes with degree 0)
    numisolnodes.g<-length(which(deg.g==0))


    ####MST-based


    # Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
    mst<-mst(g, algorithm="prim")

    # let's see if we have communities here using the
    # Grivan-Newman algorithm
    # 1st we calculate the edge betweenness, merges, etc...
    mst.communities <- edge.betweenness.community(mst, weights=NULL, directed=FALSE)
    mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
    V(mst)$color <- mst.communities$membership

    # Check the vertex degree, i.e., number of connections to each vertex
    deg.mst<-degree(mst)

    ##number of edges
    edgenum.mst<-gsize(mst)

    ###number of clusters
    clusnum.mst<-length(unique(mst.communities$membership))

    ##number of isolated nodes (i.e. nodes with degree 0)
    numisolnodes.mst<-length(which(deg.mst==0))

    output<-list(cm,g,g.communities,g.clustering,deg.g,edgenum.g,clusnum.g,numisolnodes.g,mst,mst.communities,mst.clustering,deg.mst,edgenum.mst,clusnum.mst,numisolnodes.mst)
    names(output)<-c("adjacency matrix","graph","graph communities","graph clustering","graph vertex degrees","graph number of edges","graph number of clusters","Graph number of isolated nodes","MST","MST communities","MST clustering","MST vertex degrees","MST number of edges","MST number of clusters","MST number of isolated nodes")
    return(output)

  }
}


###2.) Network differneces

######(i) overall network difference characteristics

#' Frobenius metric
#' @description The method calculates the Frobenius metric.
#' @param A,B are matrices.
#' @details
#' One of the methods to compare two networks on their overall network characteristics.
#' @return The frobenius metric of the two adjacency matrices \code{A} and \code{B}.
#' @examples
#' frobenius.metric(matrix(3,2), matrix(2,2))
#' @export
#'
frobenius.metric<-function(A,B){
  stopifnot("A and B need to be matrices or tibbles" = (class(A)==c("matrix","array")||class(A)==c("tbl_df","tbl","data.frame")) &&
                                                       (class(B)==c("matrix","array")||class(B)==c("tbl_df","tbl","data.frame")),
            "A and B need to have the same dimensions" = nrow(A)==nrow(B) && ncol(A)==ncol(B))

  output<-sqrt(sum(abs((A-B)^2)))
  return(output)
}

#' Maximum metric
#' @description The method calculates the Maximum metric.
#' @param A,B matices.
#' @details
#' One of the methods to compare two networks on their overall network characteristics.
#' @return The Maximum metric of the two adjacency matrices \code{A} and \code{B}.
#' @examples
#' max.metr(matrix(3,2), matrix(2,2))
#' @export
#'
max.metr<-function(A,B){
  stopifnot("A and B need to be matrices or tibbles" = (class(A)==c("matrix","array")||class(A)==c("tbl_df","tbl","data.frame")) &&
                                                       (class(B)==c("matrix","array")||class(B)==c("tbl_df","tbl","data.frame")),
            "A and B need to have the same dimensions" = nrow(A)==nrow(B) && ncol(A)==ncol(B))

  output<-max(abs(A-B))
  return(output)
}


#' Spectral distance
#' @description The method calculates the Spectral distance.
#' @param A,B are matrices.
#' @details
#' One of the methods to compare two networks on their overall network characteristics.
#' @return The spectral distance of the two adjacency matrices \code{A} and \code{B}.
#' @examples
#' spec.dist(matrix(3,2), matrix(2,2))
#' @export
#'
spec.dist<-function(A,B){
  stopifnot("A and B need to be matrices or tibbles" = (class(A)==c("matrix","array")||class(A)==c("tbl_df","tbl","data.frame")) &&
                                                       (class(B)==c("matrix","array")||class(B)==c("tbl_df","tbl","data.frame")),
            "A and B need to have the same dimensions" = nrow(A)==nrow(B) && ncol(A)==ncol(B))

  output<-sqrt(sum(((eigen(A)$values)-(eigen(B)$values))^2))
  return(output)
}


#' Global strength
#' @description The method calculates the Global strength.
#' @param A,B are matrices.
#' @details
#' One of the methods to compare two networks on their overall network characteristics.
#' @return The global strength of the two adjacency matrices \code{A} and \code{B}.
#' @examples
#' global.str(matrix(3,2), matrix(2,2))
#' @export
#'
global.str<-function(A,B){
  stopifnot("A and B need to be matrices or tibbles" = (class(A)==c("matrix","array")||class(A)==c("tbl_df","tbl","data.frame")) &&
                                                       (class(B)==c("matrix","array")||class(B)==c("tbl_df","tbl","data.frame")),
            "A and B need to have the same dimensions" = nrow(A)==nrow(B) && ncol(A)==ncol(B))

  output<-abs(sum(abs(A))-sum(abs(B)))
  return(output)
}


##########################characteristics used by Asada et al. (2016)

#' Differneces between two graphs
#' @description The method calculates differences between two graphs.
#' @param A,B are results from the function create.graph.
#' @details
#' One of the methods to compare two networks on their overall network characteristics.
#' @return A list of the differences in the number of edges, clusters and isolated nodes of the graph and the corresponding
#' minimum spanning tree.
#' @examples
#' diff.num(create.graph(A, list("Spearman")),create.graph(B, list("Spearman")))
#' @export
#'
diff.num<-function(A,B){

  stopifnot("A and B need to be adjacency matrices." = class(A[[1]])==c("matrix","array") && class(B[[1]])==c("matrix","array")
                                                       && nrow(A[[1]])==ncol(A[[1]]) && nrow(B[[1]])==ncol(B[[1]]),
            "The listelements 6, 7, 8, 13, 14 and 15 need to be numbers. (This is the case if A and B are results of the function create.graph.)" =
                  (class(A[[6]])=="numeric" || class(A[[6]])=="integer") && (class(A[[7]])=="numeric" || class(A[[7]])=="integer") &&
                  (class(A[[8]])=="numeric" || class(A[[8]])=="integer") && (class(B[[6]])=="numeric" || class(B[[6]])=="integer") &&
                  (class(B[[7]])=="numeric" || class(B[[7]])=="integer") && (class(B[[8]])=="numeric" || class(B[[8]])=="integer") &&
                  (class(A[[13]])=="numeric" || class(A[[13]])=="integer") && (class(A[[14]])=="numeric" || class(A[[14]])=="integer") &&
                  (class(A[[15]])=="numeric" || class(A[[15]])=="integer") && (class(B[[13]])=="numeric" || class(B[[13]])=="integer") &&
                  (class(B[[14]])=="numeric" || class(B[[14]])=="integer") && (class(B[[15]])=="numeric" || class(B[[15]])=="integer"))

  if(all(A[[1]]==0) | all(B[[1]]==0)){
    output<-rep(NA,6)
    names(output)<-c("Graph Diff num of edges","Graph Diff num of clusters","Graph diff num of isolated nodes","MST Diff num of edges","MST Diff num of clusters","MST diff num of isolated nodes")
    return(output)
  }else{

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
}


######(ii) node-specifc network difference characteristics

#' Degree differneces
#' @description The method calculates the differences between the degrees of the nodes of two matrices.
#' @param X,Y are matrices
#' @details
#' One of the methods to compare two networks on their node-specific differences.
#' @return A vector of the differneces of the degrees of the same nodes in two different networks.
#' @examples
#' degree.inv(X,Y)
#' @export
#'
degree.inv<-function(X,Y){

  stopifnot("X and Y need to be matrices" = class(X)==c("matrix","array") && class(Y)==c("matrix","array"),
            "X and Y need to have the same dimensions" = nrow(X)==nrow(Y) && ncol(X)==ncol(Y),
            "To compare the graphs of X and Y correctly they need the same columntitles" = colnames(X)==colnames(Y))


  if(all(X==0)){
    A<-graph.adjacency(X,weighted=TRUE, mode="undirected", diag=FALSE)
  }else{
    A<-graph.adjacency(X,weighted=TRUE, mode="undirected", diag=FALSE)
    # Simplfy the adjacency object
    A<-simplify(A, remove.multiple=TRUE, remove.loops=TRUE)
    # Convert edge weights to absolute values
    E(A)$weight <- abs(E(A)$weight)
  }

  if(all(Y==0)){
    B<-graph.adjacency(Y,weighted=TRUE, mode="undirected", diag=FALSE)
  }else{
    B<-graph.adjacency(Y,weighted=TRUE, mode="undirected", diag=FALSE)
    # Simplfy the adjacency object
    B<-simplify(B, remove.multiple=TRUE, remove.loops=TRUE)
    # Convert edge weights to absolute values
    E(B)$weight <- abs(E(B)$weight)
  }

  degA<-degree(A)
  degB<-degree(B)

  output<-abs(degA-degB)
  names(output)<-colnames(X)

  return(output)
}


#' Betweenness
#' @description A
#' @param X
#' @param Y
#' @details
#' One of the methods to compare two networks on their node-specific differences.
#' @return Return
#' @examples
#' betweenness.inv(X,Y)
#' @export
#'
betweenness.inv<-function(X,Y){

  stopifnot("X and Y need to be matrices" = class(X)==c("matrix","array") && class(Y)==c("matrix","array"),
            "X and Y need to have the same dimensions" = nrow(X)==nrow(Y) && ncol(X)==ncol(Y),
            "To compare the graphs of X and Y correctly they need the same columntitles" = colnames(X)==colnames(Y))


  if(all(X==0)){
    A<-graph.adjacency(X,weighted=TRUE, mode="undirected", diag=FALSE)
  }else{
    A<-graph.adjacency(X,weighted=TRUE, mode="undirected", diag=FALSE)
    # Simplfy the adjacency object
    A<-simplify(A, remove.multiple=TRUE, remove.loops=TRUE)
    # Convert edge weights to absolute values
    E(A)$weight <- abs(E(A)$weight)
  }

  if(all(Y==0)){
    B<-graph.adjacency(Y,weighted=TRUE, mode="undirected", diag=FALSE)
  }else{
    B<-graph.adjacency(Y,weighted=TRUE, mode="undirected", diag=FALSE)
    # Simplfy the adjacency object
    B<-simplify(B, remove.multiple=TRUE, remove.loops=TRUE)
    # Convert edge weights to absolute values
    E(B)$weight <- abs(E(B)$weight)
  }

  betA<-betweenness(A)
  betB<-betweenness(B)

  output<-abs(betA-betB)
  names(output)<-colnames(X)

  return(output)
}


#' Closeness
#' @description A
#' @param X
#' @param Y
#' @details
#' One of the methods to compare two networks on their node-specific differences.
#' @return Return
#' @examples
#' closeness.inv(X,Y)
#' @export
#'
closeness.inv<-function(X,Y){

  stopifnot("X and Y need to be matrices" = class(X)==c("matrix","array") && class(Y)==c("matrix","array"),
            "X and Y need to have the same dimensions" = nrow(X)==nrow(Y) && ncol(X)==ncol(Y),
            "To compare the graphs of X and Y correctly they need the same columntitles" = colnames(X)==colnames(Y))

  if(all(X==0)){
    A<-graph.adjacency(X,weighted=TRUE, mode="undirected", diag=FALSE)
  }else{
    A<-graph.adjacency(X,weighted=TRUE, mode="undirected", diag=FALSE)
    # Simplfy the adjacency object
    A<-simplify(A, remove.multiple=TRUE, remove.loops=TRUE)
    # Convert edge weights to absolute values
    E(A)$weight <- abs(E(A)$weight)
  }

  if(all(Y==0)){
    B<-graph.adjacency(Y,weighted=TRUE, mode="undirected", diag=FALSE)
  }else{
    B<-graph.adjacency(Y,weighted=TRUE, mode="undirected", diag=FALSE)
    # Simplfy the adjacency object
    B<-simplify(B, remove.multiple=TRUE, remove.loops=TRUE)
    # Convert edge weights to absolute values
    E(B)$weight <- abs(E(B)$weight)
  }

  cloA<-closeness(A)
  cloB<-closeness(B)

  output<-abs(cloA-cloB)
  names(output)<-colnames(X)

  return(output)
}


#' Eigen
#' @description A
#' @param X
#' @param Y
#' @details
#' One of the methods to compare two networks on their node-specific differences.
#' @return Return
#' @examples
#' eigen.inv(X,Y)
#' @export
#'
eigen.inv<-function(X,Y){

  stopifnot("X and Y need to be matrices" = class(X)==c("matrix","array") && class(Y)==c("matrix","array"),
            "X and Y need to have the same dimensions" = nrow(X)==nrow(Y) && ncol(X)==ncol(Y),
            "To compare the graphs of X and Y correctly they need the same columntitles" = colnames(X)==colnames(Y))

  if(all(X==0)){
    A<-graph.adjacency(X,weighted=TRUE, mode="undirected", diag=FALSE)
  }else{
    A<-graph.adjacency(X,weighted=TRUE, mode="undirected", diag=FALSE)
    # Simplfy the adjacency object
    A<-simplify(A, remove.multiple=TRUE, remove.loops=TRUE)
    # Convert edge weights to absolute values
    E(A)$weight <- abs(E(A)$weight)
  }

  if(all(Y==0)){
    B<-graph.adjacency(Y,weighted=TRUE, mode="undirected", diag=FALSE)
  }else{
    B<-graph.adjacency(Y,weighted=TRUE, mode="undirected", diag=FALSE)
    # Simplfy the adjacency object
    B<-simplify(B, remove.multiple=TRUE, remove.loops=TRUE)
    # Convert edge weights to absolute values
    E(B)$weight <- abs(E(B)$weight)
  }

  ecA<-eigen_centrality(A)$vector
  ecB<-eigen_centrality(B)$vector

  output<-abs(ecA-ecB)
  names(output)<-colnames(X)

  return(output)
}


######(iii) edge-specific network difference characteristics

##a. "considers directions (signs) of correlations"

#' Edge with direction
#' @description A
#' @param A,B are matrices.
#' @details
#' One of the methods to compare two networks on their edge-specific differneces.
#' @return Return
#' @examples
#' edge.inv.direc(A,B)
#' @export
#'
edge.inv.direc<-function(A,B){

  stopifnot("A and B need to be matrices or tibbles" = (class(A)==c("matrix","array")||class(A)==c("tbl_df","tbl","data.frame")) &&
                                                       (class(B)==c("matrix","array")||class(B)==c("tbl_df","tbl","data.frame")),
            "A and B need to have the same dimensions" = nrow(A)==nrow(B) && ncol(A)==ncol(B))

  output<-abs(A-B)
  rownames(output)<-rownames(A)
  colnames(output)<-colnames(A)
  return(output)
}


##b. "does not consider directions (signs) of correlations"

#' Edge without direction
#' @description A
#' @param A a matrix.
#' @param B a matrix.
#' @details
#' One of the methods to compare two networks on their edge-specific differneces.
#' @return Return
#' @examples
#' edge.inv(A,B)
#' @export
#'
edge.inv<-function(A,B){

  stopifnot("A and B need to be matrices or tibbles" = (class(A)==c("matrix","array")||class(A)==c("tbl_df","tbl","data.frame")) &&
                                                       (class(B)==c("matrix","array")||class(B)==c("tbl_df","tbl","data.frame")),
            "A and B need to have the same dimensions" = nrow(A)==nrow(B) && ncol(A)==ncol(B))

  output<-abs(abs(A)-abs(B))
  rownames(output)<-rownames(A)
  colnames(output)<-colnames(A)
  return(output)
}

###3.) Permutation test

#' Permutation test network
#' @description A
#' @param A a matrix.
#' @param B a matrix.
#' @param permnum the number of permutations.
#' @param methodlist the methods used to create the adjacency matrices.
#' @param thresh
#' @param score.funct the function used to compare the adjacency matrices A and B.
#' @param paired
#' @details
#' Details
#' @return Return
#' @examples
#' perm.test.nw(A = x1, B = x2, permnum = 10, methodlist = list("DistCorr"), thresh = 0.05, score.funct = frobenius.metric)
#' @export
#'
perm.test.nw <- function(A, B, permnum, methodlist, thresh, score.funct, paired=FALSE){

  stopifnot("A and B need to be matrices or tibbles with the same dimensions." = (class(A)==c("matrix","array")||class(A)==c("tbl_df","tbl","data.frame")) &&
                                                                                 (class(B)==c("matrix","array")||class(B)==c("tbl_df","tbl","data.frame")) &&
                                                                                 nrow(A)==nrow(B) && ncol(A)==ncol(B),
            "permnum needs to be a natural number." = permnum > 0 && permnum%%1 == 0,
            "thresh needs to be a non-negative number." = thresh >= 0,
            "score.funct needs to be a function." = class(score.funct)=="function",
            "paired needs to be a boolean." = paired==TRUE || paired==FALSE)

  ###compute score for original unpermuted raw data
  if(identical(score.funct, diff.num)){
    ###compute score for original unpermuted raw data
    graph.A<-create.graph(A,methodlist=methodlist,thresh=thresh)
    graph.B<-create.graph(B,methodlist=methodlist,thresh=thresh)
    value.orig<-diff.num(graph.A,graph.B)
    names(value.orig)<-c("Graph Diff num of edges","Graph Diff num of clusters","Graph Diff num of isolated nodes","MST Diff num of edges","MST Diff num of clusters","MST Diff num of isolated nodes")

  }else{
    adj.A<-create.adjacency.matrix(A,methodlist=methodlist,thresh=thresh)
    adj.B<-create.adjacency.matrix(B,methodlist=methodlist,thresh=thresh)
    value.orig<-score.funct(adj.A,adj.B)
  }

  if(paired==FALSE){

    D<-rbind(A,B)

    #####function to create permutations from raw data
    create.permut.unpaired<-function(seedex){
      set.seed(seedex)
      s<-sample(nrow(D),replace=FALSE)
      D.new<-D[s,]
      A.new<-D.new[1:nrow(A),]
      B.new<-D.new[(nrow(A)+1):(nrow(D.new)),]
      output<-list(A.new,B.new)
      return(output)
    }

    ##create permutations for raw data
    shuffle <- lapply(1:permnum, create.permut.unpaired)

  }


  if(paired==TRUE){

    #####function to create permutations from raw data
    create.permut.paired<-function(seedex){
      set.seed(seedex)
      s<-sample(c("A","B"),nrow(A),replace=TRUE)
      A.new<-rbind(A[s=="A",],B[s=="B",])
      B.new<-rbind(B[s=="A",],A[s=="B",])
      output<-list(A.new,B.new)
      return(output)
    }

    ##create permutations for raw data
    shuffle <- lapply(1:permnum, create.permut.paired)

  }


  if(any(identical(score.funct, frobenius.metric),identical(score.funct, max.metr),identical(score.funct, spec.dist),identical(score.funct, global.str))){
    ##compute score for the permnum permuted data sets
    value.perm<-rep(NA,permnum)
    for (n in 1:permnum){
      adj.A.perm<-create.adjacency.matrix(shuffle[[n]][[1]],methodlist=methodlist,thresh=thresh)
      adj.B.perm<-create.adjacency.matrix(shuffle[[n]][[2]],methodlist=methodlist,thresh=thresh)
      value.perm[n]<-score.funct(adj.A.perm,adj.B.perm)
    }

    ###permutation pseudocount p-value
    num.extr<-sum(value.perm>=value.orig)
    pvalue.ecdf<-num.extr/permnum
    pvalue.ecdf.pseudo<-(1+num.extr)/(permnum+1)

    output<-list(adj.A,adj.B,value.orig,value.perm,pvalue.ecdf,pvalue.ecdf.pseudo)
    names(output)<-c("adjancency.matrix.A","adjacency.matrix.B","test.statistic","test.statistics.perm","pvalue","pvalue.pseudocount")

  }


  if(identical(score.funct,diff.num)){
    ##compute score for the permnum permuted data sets
    value.perm<-array(data=NA,dim=c(permnum,length(value.orig)))
    for (n in 1:permnum){
      graph.A.perm<-create.graph(shuffle[[n]][[1]],methodlist=methodlist,thresh=thresh)
      graph.B.perm<-create.graph(shuffle[[n]][[2]],methodlist=methodlist,thresh=thresh)
      value.perm[n,]<-diff.num(graph.A.perm,graph.B.perm)
    }

    ###permutation pseudocount p-value
    pvalue.ecdf<-rep(NA,length(value.orig))
    pvalue.ecdf.pseudo<-rep(NA,length(value.orig))
    for (j in 1:length(value.orig)){

      if(all(is.na(value.orig[j])==TRUE)){
        pvalue.ecdf[j]<-NA
        pvalue.ecdf.pseudo[j]<-NA
      }

      else{

        if (all(is.na(value.perm[,j])==FALSE)){
          num.extr<-sum(value.perm[,j]>=value.orig[j])
          pvalue.ecdf[j]<-num.extr/permnum
          pvalue.ecdf.pseudo[j]<-(1+num.extr)/(permnum+1)
        }else{
          permnum<-length(na.omit(value.perm[,j]))
          num.extr<-sum(na.omit(value.perm[,j])>=value.orig[j])
          pvalue.ecdf[j]<-num.extr/permnum
          pvalue.ecdf.pseudo[j]<-(1+num.extr)/(permnum+1)
        }

      }

    }
    names(pvalue.ecdf)<-c("Graph Diff num of edges","Graph Diff num of clusters","Graph Diff num of isolated nodes","MST Diff num of edges","MST Diff num of clusters","MST Diff num of isolated nodes")
    names(pvalue.ecdf.pseudo)<-c("Graph Diff num of edges","Graph Diff num of clusters","Graph Diff num of isolated nodes","MST Diff num of edges","MST Diff num of clusters","MST Diff num of isolated nodes")


    output<-list(graph.A,graph.B,value.orig,value.perm,pvalue.ecdf,pvalue.ecdf.pseudo)
    names(output)<-c("graph.A","graph.B","test.statistic","test.statistics.perm","pvalue","pvalue.pseudocount")

  }


  if(any(identical(score.funct,degree.inv),identical(score.funct,betweenness.inv),identical(score.funct,closeness.inv),identical(score.funct,eigen.inv))){
    ##compute score for the permnum permuted data sets
    value.perm<-array(data=NA,dim=c(permnum,length(value.orig)))
    for (n in 1:permnum){
      adj.A.perm<-create.adjacency.matrix(shuffle[[n]][[1]],methodlist=methodlist,thresh=thresh)
      adj.B.perm<-create.adjacency.matrix(shuffle[[n]][[2]],methodlist=methodlist,thresh=thresh)
      value.perm[n,]<-score.funct(adj.A.perm,adj.B.perm)
    }

    ###permutation pseudocount p-value
    pvalue.ecdf<-rep(NA,length(value.orig))
    pvalue.ecdf.pseudo<-rep(NA,length(value.orig))

    names(pvalue.ecdf)<-names(value.orig)
    names(pvalue.ecdf.pseudo)<-names(value.orig)

    for(i in 1:length(value.orig)){
      num.extr<-sum(value.perm[,i]>=value.orig[i])
      pvalue.ecdf[i]<-num.extr/permnum
      pvalue.ecdf.pseudo[i]<-(1+num.extr)/(permnum+1)
    }

    output<-list(adj.A,adj.B,value.orig,value.perm,pvalue.ecdf,pvalue.ecdf.pseudo) #,p.ecdf,p.ecdf.pseudo)
    names(output)<-c("adjancency.matrix.A","adjacency.matrix.B","test.statistic","test.statistics.perm","pvalue","pvalue.pseudocount")

  }


  if(any(identical(score.funct, edge.inv),identical(score.funct,edge.inv.direc))){
    ##compute score for the permnum permuted data sets
    value.perm<-array(data=NA,dim=c(permnum,dim(adj.A)))
    for (n in 1:permnum){
      adj.A.perm<-create.adjacency.matrix(shuffle[[n]][[1]],methodlist=methodlist,thresh=thresh)
      adj.B.perm<-create.adjacency.matrix(shuffle[[n]][[2]],methodlist=methodlist,thresh=thresh)
      value.perm[n,,]<-score.funct(adj.A.perm,adj.B.perm)
    }

    ###permutation pseudocount p-value
    pvalue.ecdf<-array(data=NA,dim=dim(value.orig),dimnames=list(rownames(value.orig),colnames(value.orig)))
    pvalue.ecdf.pseudo<-array(data=NA,dim=dim(value.orig),dimnames=list(rownames(value.orig),colnames(value.orig)))

    for(i in 1:dim(value.orig)[1]){
      for(j in 1:dim(value.orig)[2]){
        num.extr<-sum(value.perm[,i,j]>=value.orig[i,j])
        pvalue.ecdf[i,j]<-num.extr/permnum
        pvalue.ecdf.pseudo[i,j]<-(1+num.extr)/(permnum+1)
      }
    }

    p.ecdf<-data.frame(var1=rownames(pvalue.ecdf)[row(pvalue.ecdf)[upper.tri(pvalue.ecdf)]],
                       var2=colnames(pvalue.ecdf)[col(pvalue.ecdf)[upper.tri(pvalue.ecdf)]],
                       pvalue=pvalue.ecdf[upper.tri(pvalue.ecdf)])


    p.ecdf.pseudo<-data.frame(var1=rownames(pvalue.ecdf.pseudo)[row(pvalue.ecdf.pseudo)[upper.tri(pvalue.ecdf.pseudo)]],
                              var2=colnames(pvalue.ecdf.pseudo)[col(pvalue.ecdf.pseudo)[upper.tri(pvalue.ecdf.pseudo)]],
                              pvalue=pvalue.ecdf.pseudo[upper.tri(pvalue.ecdf.pseudo)])


    output<-list(adj.A,adj.B,value.orig,value.perm,pvalue.ecdf,pvalue.ecdf.pseudo,p.ecdf,p.ecdf.pseudo)
    names(output)<-c("adjancency.matrix.A","adjacency.matrix.B","test.statistic","test.statistics.perm","pvalue","pvalue.pseudocount","pvalue.summ","pvalue.pseudocount.summ")

  }


  return(output)

}
