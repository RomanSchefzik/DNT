# DNT

Differential network testing

## Installation

### Prerequisites

Please make sure that you have the latest version of R installed.  

### Installation of the DNT package

To install the DNT package, run the following:  
`install.packages("remotes")`  
`remotes::install_github("RomanSchefzik/DNT")` 

The DNT package can then be loaded using  
`library(DNT)`

## Usage

The package provides functions to test for differences between two statistical networks based on permutation tests. There are various options for network estimation and network comparison (using various types of overall, node-specifc or edge-specific network difference chacteristics). Moreover, tools for graphical illustrations and comparisons are offered.

The core function of the package is `perm.nw.test`, which allows to test for differences between two statistical networks, where several settings have to be specified. 
