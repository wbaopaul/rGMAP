# GMAP: identifying TADs and subTADs from Hi-C data

GMAP is an algorithm to call topologically associating domains (TAD) and subdomains (subTAD) from normalized Hi-C data.
It's implemented through a R package rGMAP.


## Install from Github 
```
library(devtools)
install_github("wbaopaul/rGMAP")
```

## Install from source codes

* Download source codes [here](https://www.dropbox.com/s/f2p9x7r98u6d285/rGMAP_1.4.tar.gz?dl=0) 
and In R type:
 
```
install.packages('path to rGMAP_1.4.tar.gz', type = 'source', rep = NULL)
```

## Note
* The latest *rGMAP* was built under R-3.5.1

## Input
* For a single chromosome, a HiC contact matrix *hic_obj* supports three types of format: 
  1. a 3-column Hi-C contact matrix, corresponding to the i_th, j_th bin of a chromosom and the contact number; 
  2. a n by n matrix, with (i,j) th element corresponding to contact number between the i_th and j_th bin of a chromosome;
  3. a tab or space delimited text file of the above two types of data
  
* For multiple chromosomes, a genomic coordinate index file *index_obj* for HiC bin was required, and *hic_obj* and *index_obj* are compatible with HiC-Pro stype HiC matrix and index files. Both *hic_obj* and *index_obj* supports R data frame/data table and tab/space delimited files

  - An example of *index_obj (chromosome start end id)* in 10kb resolution:

  ```
  chr1	0	10000	1
  chr1	10000	20000	2
  chr1	20000	30000	3
  ......
  ```

  - An example of corresponding 3-column *hic_mat* file (*bin_i bin_j count*):

  ```
  10	11	1.15
  10	15	1.89
  15	20	2.20
  ......
  ```


## Output
  * data frames providing the genomic coordinates of the identified hierarchical domains
  * the final parameters for calling TADs

## Vignette
* Detailed [vignette](https://www.dropbox.com/s/n0bsr80fvmi1tp4/rGMAP-vignette.html?dl=0) for the latest version 1.4.

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
pp = plotdom(hic_rao_IMR90_chr15, NULL, res$hierTads, NULL, 5950, 6950, 30, 10000)
pp$p2



## for more information of usage of plotdom
help(plotdom)





```

## Reference
The detailed information of GMAP algorithm is described in the following paper:

[Yu, W., He, B., & Tan, K. (2017). Identifying topologically associating domains and subdomains by Gaussian Mixture model And Proportion test. Nature Communications, 8, 535. ](http://doi.org/10.1038/s41467-017-00478-8)


