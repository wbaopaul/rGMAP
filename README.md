# GMAP: identifying TADs and subTADs from Hi-C data

GMAP is an algorithm to call topologically associating domains (TAD) and subdomains (subTAD) from normalized Hi-C data.
It's implemented through a R package rGMAP.


## Install from Github 
```
library(devtools)
install_github("wbaopaul/rGMAP")
```

## Install from source codes

Download source codes [here](https://www.dropbox.com/sh/27es1vimtf5745t/AADLhBXE_wgrUIlnDS0LWpqYa?dl=0) (previous versions are also included)
```
install.packages('path to rGMAP_1.3.tar.gz', type = 'source', rep = NULL)
```
## Usage
* Input:
  - The Input is either a 3 columns (<i> <j> <count>) Hi-C map for a given chromosome, corrsponding to the i_th bin, j_th bin for the chromosome and the contact number for a contact; **note that the first two columns should be integer and in units of bin.**
  - Or a n by n contact matrix, n is the total number of bins for a chromosome

* Output:
  - data frames providing the coordinates of the identified hierarchical domains
  - the final parameters for calling TADs

## Vignette
* Detailed [vignette](https://www.dropbox.com/s/n0bsr80fvmi1tp4/rGMAP-vignette.html?dl=0) for the latest version 1.3.

## A quick example
* A quick instruction and example:

```
library(rGAMP)
help(rGAMP)

## use an example data from Rao et al. (2014 Cell)
hic_rao_IMR90_chr15   # normalized Hi-C data for IMR90, chr15 with resolution 10kb
res = rGMAP(hic_rao_IMR90_chr15, resl = 10 * 1000)
names(res)


## quickly visualize some hierarchical domains
pp = plotdom(hic_rao_IMR90_chr15, res$hierTads, 5950, 6950, 30, 10000)
pp$p2



## for more information of usage of plotdom
help(plotdom)





```

## Reference
The detailed information of GMAP algorithm is described in the following paper:

[Yu, W., He, B., & Tan, K. (2017). Identifying topologically associating domains and subdomains by Gaussian Mixture model And Proportion test. Nature Communications, 8, 535. ](http://doi.org/10.1038/s41467-017-00478-8)


