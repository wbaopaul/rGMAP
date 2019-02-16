## ----setup, include = FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = F, message = F,
  comment = "#>"
)

## ---- eval=F---------------------------------------------------------------
#  library(devtools)
#  install_github('wbaopaul/rGMAP')
#  

## ---- eval=F---------------------------------------------------------------
#  
#  install.packages('path to rGMAP_1.4.tar.gz', repos = NULL, type = 'source')
#  

## ---- fig.height=3, fig.width=15, fig.align='center'-----------------------
library(rGMAP)
hic_rao_IMR90_chr15   
res = rGMAP(hic_rao_IMR90_chr15, index_obj = NULL, resl = 10 * 1000, dom_order = 2)
names(res)

#quickly visualize some hierarchical domains
pp = plotdom(hic_rao_IMR90_chr15, index_obj = NULL, res$hierTads, chr0 = NULL, 5950, 6900, 30, 10000)
pp$p2


## ---- fig.height=3, fig.width=15, fig.align='center'-----------------------
set.seed(1)
simu_res = data_simu('poisson-dist-hier')
true_domains = simu_res$hierTads
simu_mat = simu_res$hic_mat
predicted_domains = rGMAP(simu_mat, resl = 1)$hierTads

pp = plotdom(simu_mat, NULL, predicted_domains, NULL, 1, 1000, 20, resl = 1)
pp$p2

#true_domains
#predicted_domains



## --------------------------------------------------------------------------
devtools::session_info()

