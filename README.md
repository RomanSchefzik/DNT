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

The DNT package essentially contains three main functions: `create.adjacency.matrix`, `perm.nw.test` and `graph.plot`.

With the function `create.adjacency.matrix(A,methodlist,thresh)`, an adjacency matrix for an input data table `A` can be created, where the specific method for estimating the adjacency matrix can be specified in the `methodlist` argument. Currently, the following options are implemented: Spearman correlations, partial Spearman correlations, distance correlations and EBICglasso. See `?create.adjacency.matrix` for details.

### Examples
```
create.adjacency.matrix(ExDataA,methodlist=list("Spearman"),thresh=0.05)
```

The core function of the DNT package is `perm.nw.test`, which allows to test for differences between two statistical networks, where several settings have to be specified. First, a network estimation method has to be chosen. Here, several options are provided: See `?perm.test.nw` for details.

### Examples
```
perm.nw.test(ExDataA,ExDataB,permnum=10000,methodlist=list("Spearman"),thresh=0.05,score.func=global.str,paired=TRUE)
```

The function `graph.plot` provides visual comparison between two networks. See `?graph.plot` for details.
By calling `create.app()` in the R console, an R Shiny app is loaded, which implements the function `graph.plot` in a user-friendly interface.

### Examples
```
graph.plot(ExDataA,ExDataB,methodlist=list("Spearman"))
```
