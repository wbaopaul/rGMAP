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

* Output: a data frame shows the coordinates of the identified hierarchical domains

* Detailed instruction can be found by:
```
library(rGAMP)
help(rGAMP)
```

## Reference
The detailed information of GMAP algorithm is described in the following paper:

[Yu, W., He, B., & Tan, K. (2017). Identifying topologically associating domains and subdomains by Gaussian Mixture model And Proportion test. Nature Communications, 8, 535. ](http://doi.org/10.1038/s41467-017-00478-8)


