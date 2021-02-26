# DNT

Differential network tests

## Installation

### Prerequisites

Please make sure that you have the latest version of R installed.  

### Installation of the DNT package

To install the DNT package, run the following:  
`install.packages("remotes")` 

`library(remotes)`

`remotes::install_github("RomanSchefzik/DNT")` 

The DNT package can then be loaded using  
`library(DNT)`

## Usage

The DNT package provides an overall frame for testing for differences between two statistical networks based on permutation tests. There are various options for (i) network estimation (i.e. specifying the adjacency matrices) and (ii) network comparison (using various types of overall, node-specifc or edge-specific network difference chacteristics). Moreover, tools for graphical illustrations and comparisons are offered. DNT fuses scattered functions from packages such as igraph, Hmisc, energy and qgraph in an overarching frame for differential network testing. 

The DNT package essentially contains three main functions: `create.adjacency.matrix`, `perm.test.nw` and `comp.plot`.

### The function `create.adjacency.matrix`

With the function `create.adjacency.matrix(A,methodlist,thresh)`, an adjacency matrix representing a network for an input data table `A` can be created, where the specific method for estimating the adjacency matrix can be specified in the `methodlist` argument. Currently, the following options are implemented: Spearman correlations, partial Spearman correlations, distance correlations (for all these options, edge weights are set to zero if the corresponding (adjusted or unadjusted) p-value is greater than the specified significance threshold `thresh`, which is typically set to 0.05) and EBICglasso. <br/>
See `?create.adjacency.matrix` for details, in particular on how to exactly set up the `methodlist` argument.

#### Examples (using the supplied example data sets `ExDataA` and `ExDataB`)
```
create.adjacency.matrix(A=ExDataA,methodlist=list("Spearman"),thresh=0.05)
create.adjacency.matrix(A=ExDataA,methodlist=list("PCSpearman.adj","BH"),thresh=0.1)
create.adjacency.matrix(A=ExDataB,methodlist=list("EBICglasso","spearman",0.1))
```

### The function `perm.test.nw`

The main function of the DNT package is `perm.test.nw(A,B,permnum,methodlist,thresh,score.funct,paired)`, which allows to test for differences between two statistical networks derived from input data sets `A` and `B` in a testing frame based on a specified number of permutations `permnum`. Both a network estimation method and a network difference charcteristic, according to which the network difference is measured, has to be chosen. For network estimation and creation of adjacency matrices, the above-mentioned methods can be specified in the argument `methodlist` (and, accompanying, `thresh`). The considered network difference characteristic has to be specified in the `score.funct` argument, where several options for overall (Frobenius metric, difference in global strength, maximum metric, spectral distance, difference in the number of edges/clusters/isolated nodes), node-specific (difference in betweenness/closeness/degree/eigenvector centrality) or edge-specific (difference in (absolute) edge weights) network difference chracteristics are provided. Finally, in the `paired` argument, it can be chosen whether a paired or unpaired version of the permutation test frame shall be considered, depending on the structure of the input data. The output of the `perm.test.nw` function comprises the estimated adjacency matrices, the value(s) of the test statistic(s), the values of the test statistic(s) when applying permutations, the corresponding p-values and the corresponding p-values when including a pseudocount to avoid p-values that are exactly zero in permutation-based settings.<br/>
See `?perm.test.nw` for details, in particular on how to exactly set up the `methodlist` and `score.funct` arguments, and for the precise output format of the function.

#### Examples (using the supplied example data sets `ExDataA` and `ExDataB`)
```
res1<-perm.test.nw(A=ExDataA,B=ExDataB,permnum=10000,methodlist=list("Spearman"),thresh=0.05,score.func=global.str,paired=TRUE)
res2<-perm.test.nw(A=ExDataA,B=ExDataB,permnum=10000,methodlist=list("EBICglasso","spearman",0.1), score.funct=degree.inv)
res3<-perm.test.nw(A=ExDataA,B=ExDataB,permnum=10000,methodlist=list("Spearman.adj","bonferroni"),score.funct=edge.inv)
```

### The function `comp.plot`

The function `comp.plot` provides a tool for visual comparison between two networks, where the above-mentioned methods can be used for network estimation and creating the adjacency matrices. Moreover, several options for setting up the layout and the graphical representation are available. <br/>
See `?comp.plot` for details.

#### Examples (using the supplied example data sets `ExDataA` and `ExDataB`)
```
comp.plot(A=ExDataA,B=ExDataB,methodlist=list("Spearman"))
comp.plot(A=ExDataA,B=ExDataB,methodlist=list("PCSpearman.adj","bonferroni"),layout=igraph::layout.circle,curved=FALSE)
comp.plot(A=ExDataA,B=ExDataB,methodlist=list("EBICglasso","spearman",0.1),layout=igraph::layout.fruchterman.reingold,curved=FALSE)
```

![ScreenShot](/ExampleFigures/ExPlotSpearman.png)

<img src="/ExampleFigures/ExPlotSpearman.png" width="500" height="500" class="center">

### R Shiny app for visual comparison of two networks

By calling `create.app()` in the R console, an R Shiny app in which the function `graph.plot` for visual comparison of two networks is implemented in a user-friendly interface will be opened.
```
create.app()
```
