# GMAP: identifying TADs and subTADs from Hi-C experiment

GMAP is an algorithm to call topologically associating domains (TAD) and subdomains (subTAD) from normalized Hi-C data.
It's implemented through a R package *rGMAP*.


## Installation 
```
library(devtools)
install_github("wbaopaul/rGMAP")
```

## Alternative installation if you cannot install it directly from github
Download source package here https://www.dropbox.com/s/ami6hopkl6c6wlk/rGMAP_1.2.tar.gz?dl=0 and then install it in R by:
```
install.packages('rGMAP_1.2.tar.gz', type = 'source', rep = NULL)
```

## Usage
* **Input:**
  - The Input is either a 3 columns Hi-C map for a given chromosome, with each row represents the start bin index, end bin index and the normalized counts (score) for a contact. **Note the first two columns should be in the unit of bin**.
  - Or a n by n contact matrix, n is the total number of bins for a chromosome

* **Output:**
  - data frames providing the genomic coordinates of the identified hierarchical domains
  - the final parameters for calling TADs

* **Example**: 

```
library(rGMAP)
help(rGMAP)

## use an example data from Rao et al. (2014 Cell)
hic_rao_IMR90_chr15   # normalized Hi-C data for IMR90, chr15 with resolution 10kb
head(hic_rao_IMR90_chr15)

res = rGMAP(hic_rao_IMR90_chr15, resl = 10 * 1000, dom_order = 2)
names(res)


## quickly visualize some hierarchical domains
pp = plotdom(hic_rao_IMR90_chr15, res$hierTads, 6000, 6500, 30, 10)
pp$p2



## for more information of usage of plotdom
help(plotdom)





```

## Reference
The detailed information of GMAP algorithm is described in the following paper:

[Yu, W., He, B., & Tan, K. (2017). Identifying topologically associating domains and subdomains by Gaussian Mixture model And Proportion test. Nature Communications, 8, 535. ](http://doi.org/10.1038/s41467-017-00478-8)


