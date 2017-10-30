# GMAP: identifying TADs and subTADs from Hi-C data

GMAP is an algorithm to call topologically associating domains (TAD) and subdomains (subTAD) from normalized Hi-C data.
It's implemented through a R package rGMAP.


## Installation 
```
library(devtools)
install_github("wbaopaul/rGMAP")
```
## Usage
* Input:
  - The Input is either a 3 columns Hi-C map for a given chromosome, with each row corrsponding to the start bin, end bin and the contact number for a contact
  - Or a n by n contact matrix, n is the total number of bins for a chromosome

* Output:
  - data frames providing the coordinates of the identified hierarchical domains
  - the final parameters for calling TADs

* Detailed instruction and an example can be found by:

```
library(rGAMP)
help(rGAMP)

## use an example data from Rao et al. (2014 Cell)
hic_rao_IMR90_chr15   # normalized Hi-C data for IMR90, chr15 with resolution 10kb
res = rGMAP(hic_rao_IMR90_chr15, resl = 10 * 1000, dom_order = 2, bthr = 400)
names(res)


## quickly visualize some hierarchical domains
pp = plotdom(hic_rao_IMR90, res$hierTads, 6000, 6500, 20, 10)
pp$p2



## for more information of usage of plotdom
help(plotdom)





```

## Reference
The detailed information of GMAP algorithm is described in the following paper:

[Yu, W., He, B., & Tan, K. (2017). Identifying topologically associating domains and subdomains by Gaussian Mixture model And Proportion test. Nature Communications, 8, 535. ](http://doi.org/10.1038/s41467-017-00478-8)


