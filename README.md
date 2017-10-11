# GMAP

GMAP is an algorithm to call topologically associating domains (TAD) and subdomains (subTAD) from normalized Hi-C data.
It's implemented through a R package rGMAP.


## Install

## Example
library(rGMAP)
set.seed(1)
simu_dat = data_simu('poisson-dist-hier') ## generate a synthetic data
hic_mat = simu_dat$hic_mat
true_tads = simu_dat$tads_true
hier_true = simu_dat$hierTads
true_tads
hier_true
res = rGMAP(hic_mat, resl = 1, dom_order = 2, logt = T, fcthr = 0.95)
hierTads = res$hierTads  ## predicted hierarchical domains
hierTads

## Reference
Yu, W., He, B., & Tan, K. (2017). Identifying topologically associating domains and subdomains by Gaussian Mixture model And Proportion test. Nature Communications, 8, 535

