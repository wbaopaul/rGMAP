##---------------------------------------------------------------------##
## --This program is written by Wenbao Yu for implementing GMAP -------##
## --This file contains all the source codes
##---------------------------------------------------------------------##

# call upstream window counts of class1
up_stat <- function(x, pp, d = 50){
  pairx = (x - d + 1) : x
  if(x < d) pairx = 1:x
  if(x < 5) pairx = 1 : 5   # use 4 as a buffer
  num_class1 = sum(pp[pairx, pairx])
  tnum = (length(pairx))^2
  return(c(num_class1, tnum))
}
up_stat = cmpfun(up_stat)


# call downstream window counts of class1
down_stat <- function(x, pp, d = 50){
  pairx = x:(x + d - 1)
  if(x > nrow(pp) - d + 1) pairx = (nrow(pp) - d + 1) : nrow(pp)
  if(x > nrow(pp) - 5) pairx = (x-4) : nrow(pp)

  num_class1 = sum(pp[pairx, pairx])
  tnum = (length(pairx))^2
  return(c(num_class1,  tnum))
}
down_stat = cmpfun(down_stat)


# call counts of class1 beween the upwindow and downwindow
btw_stat <- function(x, pp, d = 50){
  pairx1 = (x - d + 1) : x
  pairx2 = x : (x + d - 1)

  if(x < d) pairx1 = 1:x
  if(x < 5) pairx1 = 1 : 5

  if(x > nrow(pp) - d + 1) pairx2 = (nrow(pp) - d + 1) : nrow(pp)
  if(x > nrow(pp) - 5) pairx2 = (x-4) : nrow(pp)

  num_class1 = sum(pp[pairx1, pairx2]) + sum(pp[pairx2, pairx1])
  tnum = 2 * length(pairx1) * length(pairx2)
  return(c(num_class1,  tnum))
}
btw_stat = cmpfun(btw_stat)


## calculate within-between test stat
cal_stat <- function(pp, d = 50){
  bpos = 1:nrow(pp)
  ustat = lapply(bpos, up_stat,  pp=pp, d = d)
  dstat = lapply(bpos, down_stat, pp=pp,  d = d)
  bstat = lapply(bpos, btw_stat, pp=pp,  d = d)

  n_up = sapply(ustat, function(x) x[2])
  n_btw = sapply(bstat, function(x) x[2])
  n_down = sapply(dstat, function(x) x[2])

  prop_up = sapply(ustat, function(x) x[1]/x[2])
  prop_btw = sapply(bstat, function(x) x[1]/x[2])
  prop_down = sapply(dstat, function(x) x[1]/x[2])

  prop_up[n_up == 0] = 0L
  prop_down[n_down == 0] = 0L
  prop_btw[n_btw == 0] = 0L

  # for testing difference between up and down
  p0 = (n_down * prop_down + n_up * prop_up)/(n_up + n_down)
  serror = sqrt(p0 * (1-p0) * (1/(n_up) + 1/n_down))
  test = (prop_up - prop_down)/serror
  test[serror == 0] = 100L

  stat_up = pmax(0, test)
  stat_down = pmax(0, -test)


  # testing up+down > between
  n = length(bpos)
  p0 = (n_down * prop_down + n_up * prop_up + n_btw * prop_btw) / (n_down + n_btw + n_up)
  serror = sqrt(p0 * (1-p0) * (1/(n_up + n_down) + 1/n_btw))
  prop_within = (n_down * prop_down + n_up * prop_up)/(n_up + n_down)
  stat_wb = (prop_within - prop_btw)/serror
  stat_wb[serror == 0] = 0L


  ## proportion fold change between intra and inter domains
  fc1 = prop_up/prop_btw
  fc2 = prop_down/prop_btw
  fc1[is.infinite(fc1)] = 1000
  fc2[is.infinite(fc2)] = 1000

  fc1[is.na(fc1)] = 1
  fc2[is.na(fc2)] = 1

  pr_fc_wb = pmax(fc1, fc2)


  ## proportion fold change between up and down domains
  pr_fc_ud = pmax(prop_up/prop_down, prop_down/prop_up)

  pr_fc_ud[is.infinite(pr_fc_ud)] = 1000
  pr_fc_ud[is.na(pr_fc_ud)] = 1



  return(list('up' = round(stat_up, 2), 'down' = round(stat_down,2),
              'wb' = round(stat_wb, 2),
              'pr_fc_wb' = round(pr_fc_wb, 2), 'pr_fc_ud' = round(pr_fc_ud,2)))

}
cal_stat = cmpfun(cal_stat)



## find local peaks for a given vector
## better version
localMaxima <- function(x, pr_fc_wb, thr = 2, fcthr = 0.9, dp = 10) {
  #y = which(diff(sign(diff(x))) == -2) + 1

  # smmothing using sliding window before search local peaks -- not used
  x = runmean(x, 5)

  x = round(x, 1)
  y = extrema(x)$maxindex[, 1]
  xy = x[y]
  y = y[xy >= thr]  ## screening by t1

  pr_fc_wb = round(runmean(pr_fc_wb, 5), 2)


  xy = pr_fc_wb[y]  ## further screening some peaks with small wb-proportion fold change
  y = y[xy >= quantile(pr_fc_wb, fcthr)]

  if(length(y) == 0) return(NULL)

  # remove peaks that are close
  if(!is.null(y)) y = combPeaks(y, x, dp)


  y = y[y >= 6]
  y = y[y <= (length(x) - 6)]

  return(y)
}
localMaxima = cmpfun(localMaxima)


# define fcthr,thr locally
localMaxima1 <- function(x,  pr_fc_wb, thr = 0.75, fcthr = 0.9, dp = 10) {
  # smmothing using sliding window before search local peaks
  x = runmean(x, 5)
  x = round(x, 1)

  local_t1_quantile <- function(ly){
    ll = max(ly - 500, 1)
    rr = min(ly + 500, length(x))
    return(quantile(x[ll : rr], thr))
  }

  local_t1thr = quantile(x, thr)

  y = extrema(x)$maxindex[, 1]
  xy = x[y]

  if(length(x) >= 1000) local_t1thr = sapply(y, local_t1_quantile)


  y = y[xy >= local_t1thr]  ## screening by t1

  pr_fc_wb = round(runmean(pr_fc_wb, 5), 2)

  ## locally screening further, using fold change
  xy_fc = pr_fc_wb[y]
  local_fcthr = quantile(pr_fc_wb, fcthr)

  local_fc_quantile <- function(ly){
    ll = max(ly - 500, 1)
    rr = min(ly + 500, length(x))
    return(quantile(pr_fc_wb[ll : rr], fcthr))
  }

  if(length(x) >= 1000) local_fcthr = sapply(y, local_fc_quantile)
  y = y[xy_fc >= local_fcthr]

  if(length(y) == 0) return(NULL)

  # remove peaks that are close
  if(!is.null(y)) y = combPeaks(y, x, dp)


  y = y[y >= 6]
  y = y[y <= (length(x) - 6)]

  return(y)
}
localMaxima1 = cmpfun(localMaxima1)


## define tad by gmap given peaks and statistics
def_tad_gmap <- function(wb_peaks, stat_up, stat_down, nn, thr = 2, pr_fc_ud){
  s1 = stat_up[wb_peaks]
  s2 = stat_down[wb_peaks]
  fc = pr_fc_ud[wb_peaks]

  len = length(wb_peaks)
  signs = rep(0, len)
  pr_fc_ud = pr_fc_ud[wb_peaks]
  ind1 = which(s1 > thr & fc > quantile(pr_fc_ud, 0.95))
  ind2 = which(s2 > thr & fc > quantile(pr_fc_ud, 0.95))
  if(length(ind1) > 0) signs[ind1] = 1
  if(length(ind2) > 0) signs[ind2] = -1

  start = end = NULL


  # start call first tad
  if(signs[1] == -1){
    start = c(start, 1)
    i = 1
  }
  if(signs[1] == 0){
    start = c(start, 0)
    end = c(end, 1)
    start = c(start, 1)
    i = 1
  }
  if(signs[1] == 1){
    start = c(start, 0)
    i = 0
  }

  repeat{
    if(i == len) break
    for(j in (i+1):len){

      if(signs[j] == -1 & j != len) next
      if(signs[j] == -1 & j == len) break

      if(signs[j] == 0){
        end = c(end, j)
        start = c(start, j)
        break
      }

      if(signs[j] == 1){
        if(j==len){
          end = c(end, j)
          break
        }

        if(j < len & signs[j+1] != 1) {
          end = c(end, j)
          start = c(start, j+1)
          break
        }

        if(j < len & signs[j+1] == 1) next

      }
    }  ## end of for loop
    i = max(start)
    if(j == len) break
  }    ## end of repeat loop

  # final ends
  if(length(start) > length(end)){
    end <- c(wb_peaks[end], nn)
  }else{
    end <- wb_peaks[end]
  }

  # final starts
  if(start[1] == 0){
    start = c(1, wb_peaks[start[-1]])
  }else{
    start = wb_peaks[start]
  }

  dd = data.frame("start" = start, "end" = end)
  dd = dd[dd$start != dd$end, ]
  return(dd)
}
def_tad_gmap = cmpfun(def_tad_gmap)



## using boundary to define segmatation
full_seg <- function(bounds){
  len = length(bounds)
  return(data.frame('start' = bounds[-len], 'end' = bounds[-1]))
}
full_seg = cmpfun(full_seg)


## score function
## score = difference between TAD area and a predefined area (not whole background)
## the predifined area: wd bins from diagonal of hic_mat
tune_score <- function(pp, tads, rnum_bk, tnum_bk){

  # calculate real contact number and total number for within TADs
  rnum_tad = apply(as.matrix(tads), 1, function(x) sum(pp[x[1]:x[2], x[1]:x[2]]))
  tnum_tad = apply(as.matrix(tads), 1, function(x) (x[2]-x[1]+1)^2)

  # define the score
  rnum_tad = sum(rnum_tad)
  tnum_tad = sum(tnum_tad)

  p1 = rnum_tad/tnum_tad
  p2 = rnum_bk/tnum_bk

  p0 = (rnum_tad + rnum_bk)/(tnum_tad + tnum_bk)
  dpScore = (p1 - p2)/sqrt(p0*(1-p0)*(1/tnum_tad + 1/tnum_bk))

  if(p1 == p2) dpScore = 0



  return(dpScore)
}
tune_score = cmpfun(tune_score)



## including tune t1 and t2 -- given dp and d
tune_T1T2 <- function(pp, stats, d, dp, rnum_bk, tnum_bk, fcthr, t1thr = 0.5){

  stat_up = stats$up
  stat_down = stats$down
  stat_wb = stats$wb
  pr_fc_wb = stats$pr_fc_wb
  pr_fc_ud = stats$pr_fc_ud
  nn = nrow(pp)

  # no local peaks
  if(all(diff(stat_wb) <=0) || all(diff(stat_wb) >= 0)) return(list('score' = 0, 'tads' = NULL))

  # 75% quantile is too close to 99% quantile
  if(max(stat_wb)- quantile(stat_wb, 0.75) < 1) return(list('score' = 0, 'tads' = NULL))

  sthr = 0.05
  sthr0 = round(qnorm(sthr/nn, lower.tail = F), 1)

  # give candidates for t1
  stat_wb = runmean(stat_wb, 5)

  if(F){
    cand_t1 <- quantile(stat_wb, seq(t1thr, 1.0, length = 10))
    if(all(cand_t1 < sthr0)) return(list('score' = 0, 'tads' = NULL))

    cand_t1 = unique(pmax(sthr0, round(cand_t1, 1)))

  }
  cand_t1 = seq(t1thr, 1.0, length = 10)  # new way

  # define t2 candidates globally
  tmp = c(stat_up, stat_down)
  cand_t2 <- quantile(tmp[tmp != 0], seq(0.975, 1.0, length = 10))
  tmp = qnorm(sthr, lower.tail = F)
  cand_t2 = unique(pmax(tmp, round(cand_t2, 2)))


  len1 = length(cand_t1)
  len2 = length(cand_t2)
  score = matrix(0, len1, len2)

  ## define another stats for choosing between two close peaks
  #stat_alt = cal_stat(pp, max(d-5, floor(d/2)))$wb
  #stat_alt = cal_stat(pp, max(d-5, floor(d/2)))$pr_fc_wb

  for(i in 1:len1){
    wb_peaks = localMaxima1(stat_wb, pr_fc_wb, thr = cand_t1[i], fcthr, dp = dp)

    if(is.null(wb_peaks)) {
      score[i, ] = 0
      next
    }

    if(length(wb_peaks) == 0) {
      score[i, ] = 0
      next
    }

    for(j in 1:len2){

      tads <- def_tad_gmap(wb_peaks, stat_up, stat_down, nn, thr = cand_t2[j], pr_fc_ud)

      if(nrow(tads) == 0) {
        score[i, j] = 0
      }else{
        score[i, j] <- tune_score(pp, tads, rnum_bk, tnum_bk)
      }
    }
  }

  if(max(score, na.rm = TRUE) <= 0) return(list('score' = 0, 'tads' = NULL))

  ind = which(score == max(score, na.rm = TRUE), arr.ind=TRUE)[1, ]
  wb_peaks = localMaxima1(stat_wb, pr_fc_wb, thr = cand_t1[ind[1]], fcthr, dp = dp)

  if(is.null(wb_peaks)) return(list('score' = 0, 'tads' = NULL))

  tads <- def_tad_gmap(wb_peaks, stat_up, stat_down, nn,
                       thr = cand_t2[ind[2]], pr_fc_ud)

  return(list('tads' = tads, 'score' = max(score), 't1' = round(cand_t1[ind[1]], 3),
              't2' = round(cand_t2[ind[2]], 3) , 'stats' = stats,
              'local_peaks' = wb_peaks))
}
tune_T1T2 = cmpfun(tune_T1T2)




## tune all parameters -- give candidates of d and dp
tune_allPara <- function(pp, candD, candDp, rnum_bk, tnum_bk, fcthr, t1thr = 0.5){
  lend = length(candD)
  lendp = length(candDp)
  score = matrix(0, lend, lendp)
  nn = nrow(pp)
  sthr = 0.05

  for(i in 1:lend){
    stats = cal_stat(pp, d = candD[i])
    if(all(stats$wb <= qnorm(sthr/nn, lower.tail = F))) {
      score[i, ] = 0
      next
    }
    for(j in 1:lendp){
      score[i, j] = tune_T1T2(pp, stats, d = candD[i], dp = candDp[j],
                              rnum_bk, tnum_bk, fcthr, t1thr = t1thr)$score
    }
  }

  if(max(score, na.rm = TRUE) <= 0) return(NULL)

  ind = which(score == max(score, na.rm = TRUE), arr.ind = T)[1, ]
  d = candD[ind[1]]
  dp = candDp[ind[2]]

  stats = cal_stat(pp, d = d)
  res = tune_T1T2(pp, stats, d = d, dp = dp, rnum_bk, tnum_bk, fcthr, t1thr = t1thr)
  res$d = d
  res$dp = dp
  return(res)

}
tune_allPara = cmpfun(tune_allPara)



## remove redundant domains and add false gaps
domain_correction <- function(subTads, pp, rnum_bk, tnum_bk){

  subTads = as.matrix(subTads)
  bds = sort(unique(as.vector(subTads)))

  ## filter domains with very small contacts
  subTads = as.matrix(full_seg(bds))

  s_density = round(apply(subTads, 1, function(x) mean(pp[x[1]:x[2], x[1]:x[2]])), 2)
  tsize = round(apply(subTads, 1, function(x) x[2] - x[1]))
  nn = nrow(pp)

  calDensity_expect <- function(tsize0){
    ids = lapply(1:nn, function(x) return(cbind(x, c(max(1, x - tsize0):min(nn, x + tsize0)))))
    ids = do.call('rbind', ids)
    return(round(mean(pp[ids]), 2))
  }


  e_density = sapply(tsize, calDensity_expect)

  #tmp = cbind(subTads, s_density, e_density, e_density*1.1)
  #write.table(tmp, file = 'tmp.txt', row.names = F, col.names = F, quote = F, sep = '\t')
  #message(round(rnum_bk/tnum_bk, 2))


  ids = which(s_density <= pmax(e_density, rnum_bk/tnum_bk))
  if(length(ids) > 0) subTads = subTads[-ids, ]

  return(subTads)
}
domain_correction = cmpfun(domain_correction)


## remove redundant domains and add false gaps locally
## dw represents dw up- and down- range of the current domain
domain_correction_local <- function(subTads, pp, rnum_bk, tnum_bk, dw = 5){

  subTads = as.matrix(subTads)
  bds = sort(unique(as.vector(subTads)))

  ## filter domains with very small contacts
  subTads = as.matrix(full_seg(bds))

  s_density = round(apply(subTads, 1, function(x) mean(pp[x[1]:x[2], x[1]:x[2]])), 2)
  tsize = round(apply(subTads, 1, function(x) x[2] - x[1]))
  nn = nrow(pp)

  calDensity_expect_local <- function(tsize0, dw0 = dw){
    ids = lapply(1:nn, function(x) return(cbind(x, c(max(1, x - dw0 * tsize0):min(nn, x + dw0 * tsize0)))))
    ids = do.call('rbind', ids)
    return(round(mean(pp[ids]), 2))
  }


  e_density = sapply(tsize, calDensity_expect_local)

  #tmp = cbind(subTads, s_density, e_density, e_density*1.1)
  #write.table(tmp, file = 'tmp.txt', row.names = F, col.names = F, quote = F, sep = '\t')
  #message(round(rnum_bk/tnum_bk, 2))


  ids = which(s_density <= pmax(e_density, rnum_bk/tnum_bk))
  if(length(ids) > 0) subTads = subTads[-ids, ]

  return(subTads)
}
domain_correction_local = cmpfun(domain_correction_local)



call_domain <- function(sub_mat, Max_d, min_d, Max_dp, min_dp, Bg_d, fcthr, hthr = 0.5,
                        bthr = 200, t1thr = 0.5){

  nn = nrow(sub_mat)


  # if the tad is too small (less than 10 bins), do not call subtad
  if(nn <= 10) return(list("tads" = NULL))

  # do clustering
  if(nn <= bthr){
    temp = as.vector(sub_mat[upper.tri(sub_mat, diag=T)])
    if(length(unique(temp)) < 2) return(list("tads" = NULL))

    class1 = rep(0, length(temp))

    # deal with outliers
    a3 = quantile(temp, 0.9)

    if(a3 <= 0) return(list("tads" = NULL))

    ids = which(temp >= a3)

    if(length(ids) > 0){
      class1[ids] = 1
      temp = temp[-ids]
      class0 = class1[-ids]
    }else{
      class0 = class1
    }

    if(quantile(temp, 0.9) <= 0) return(list("tads" = NULL))

    if(length(unique(temp)) >= 2){
      model2 <- Mclust(temp, G = 2, modelNames = 'E')

      if(!is.null(model2)){
        ind = which.max(model2$parameters$mean)
        posterior = model2$z[, ind]
        class0 = ifelse(posterior > hthr, 1, 0)
        rm(posterior)
      }

      if(is.null(model2)){ ## not covergent
        message('Trying using models with different variances')
        model2 <- Mclust(temp, G = 2, modelNames = 'V')
        if(!is.null(model2)) {
          ind = which.max(model2$parameters$mean)
          posterior = model2$z[, ind]
          class0 = ifelse(posterior > hthr, 1, 0)
          rm( posterior)
        }else{
          message('No models can be fitted, using 75% quantile as cutoff')
          class0 = ifelse(temp > quantile(temp, 0.75), 1, 0)
        }
      }


      # just in case:
      if(length(unique(class0)) == 1 )  class0 = ifelse(temp > quantile(temp, 0.75), 1, 0)
      rm(temp, model2)

    }


    pp = matrix(0L, nn, nn)

    if(length(ids) > 0) {
      class1[-ids] = class0
      pp[upper.tri(pp, diag=T)] = class1
    }else{
      pp[upper.tri(pp, diag=T)] = class0
    }

    pp = pp + t(pp) - diag(diag(pp))
    rm(class0, class1)
  }

  if(nn > bthr){
    ## if the hic-matrix is too large, we assume two loci with distance larger than
    ## say 1000 bins have no interactions. Data from this part will not be used for
    ## constructing the mixture models
    if(all(class(sub_mat) != "data.frame")) {
      mat = data.frame(sub_mat)
    }
    mat = as.list(mat)  ## for faster subsetting

    colIDs =  sapply(1:nn, function(x) return(c(x : min(nn, x + bthr))))
    len = sapply(colIDs, length)

    getData <- function(x){
      return(mat[[x]][colIDs[[x]]])
    }
    temp = unlist(sapply(1:nn, getData))
    rm(mat)



    rowIDs = rep(1:nn, len)
    colIDs = unlist(colIDs)

    partialIDs = cbind(rowIDs, colIDs)
    rm(rowIDs, colIDs)

    # deal with outliers
    class1 = rep(0, length(temp))

    # deal with outliers

    a3 = (quantile(temp, 0.975) )
    a1 = max(quantile(temp, 0.025), 0)

    id1 = which(temp <= a1)
    id3 = which(temp >= a3)
    ids = NULL

    if(length(id3) > 0){
      class1[id3] = 1
      ids = c(id1, id3)
      temp = temp[-ids]
      class0 = class1[-ids]
    }else{
      if(length(id1) > 0){
        ids = id1
        temp = temp[-ids]
        class0 = class1[-ids]
      }else{
        class0 = class1
      }
    }


    if(length(unique(temp)) >= 2){
      model2 <- Mclust(temp, G = 2, modelNames = 'E')

      if(!is.null(model2)){
        ind = which.max(model2$parameters$mean)
        posterior = model2$z[, ind]
        class0 = ifelse(posterior > hthr, 1, 0)
        rm(posterior)
      }

      if(is.null(model2)){ ## not covergent
        message('Trying using models with different variances')
        model2 <- Mclust(temp, G = 2, modelNames = 'V')
        if(!is.null(model2)) {
          ind = which.max(model2$parameters$mean)
          posterior = model2$z[, ind]
          class0 = ifelse(posterior > hthr, 1, 0)
          rm(posterior)
        }else{
          message('No models can be fitted, using 75% quantile as cutoff')
          class0 = ifelse(temp > quantile(temp, 0.75), 1, 0)
        }
      }


      # just in case:
      if(length(unique(class0)) == 1 )  class0 = ifelse(temp > quantile(temp, 0.75), 1, 0)
      rm(temp, model2)
    }


    pp = matrix(0L, nn, nn)

    if(length(ids) > 0) {
      class1[-ids] = class0
      pp[partialIDs] = class1
    }else{
      pp[partialIDs] = class0
    }

    pp = pp + t(pp) - diag(diag(pp))
    rm(class0, partialIDs)
  }

  ##
  diag(pp) <- 0

  # define candidates for d and dp
  #candD = seq(min_d, Max_d, by = min_d)
  candD = seq(min_d, Max_d, by = min_d)

  candDp = seq(min_dp, Max_dp, by = min_dp)

  # define the background
  if(Bg_d >= nn){

    rnum_bk <- sum(pp)
    tnum_bk <- length(pp)

    rnum_bk0 = sum(sub_mat)

  }else{
    if(Bg_d > bthr) Bg_d = bthr
    idBg = sapply(1:nn, function(x) return(cbind(x, c(max(1, x - Bg_d):min(nn, x + Bg_d)))))
    idBg = do.call('rbind', idBg)
    temp = pp[idBg]
    rnum_bk <- sum(temp)
    tnum_bk <- length(temp)

    temp0 = sub_mat[idBg]
    rnum_bk0 = sum(temp0)
  }

  if(rnum_bk >= 0.95 * tnum_bk || rnum_bk < 1/sqrt(tnum_bk) * tnum_bk) return(list("tads" = NULL))

  res = tune_allPara(pp, candD, candDp, rnum_bk, tnum_bk, fcthr, t1thr)
  subTads = res$tads


  if(!is.null(subTads)) res$tads = domain_correction_local(subTads, pp, rnum_bk, tnum_bk)
  #if(!is.null(subTads)) res$tads = domain_correction(subTads, sub_mat, rnum_bk0, tnum_bk)

  return(res)
}
call_domain = cmpfun(call_domain)





#' Detect hierarchical choromotin domains by GMAP
#' @param  hic_mat n by n matrix, n is the number of bins. Or hic_mat could by 3 columns
#' matrix or data.frame with columns: bin1, bin2, counts, in which bin1 and bin2, from 1 to m, are the bin
#' index of a hic contact
#' @param  resl The resolution (bin size), default 10kb
#' @param  min_d The minimum d (d: window size), default 25
#' @param  Max_d The maximum d (d: window size), default 2mb (or 200 bins)
#' @param min_dp The minmum dp (dp: lower bound of tad size), defalt 5
#' @param Max_dp The maximum dp (dp: lower bound of tad size), defalt 10
#' @param Bg_d The size of the background domain, dalfault 2mb or 200bins
#' @param hthr The lower bound cutoff for posterior probability, default 0.5 for TADs,
#'        0.9 for subTADs
#' @param bthr Threshold, above which a locus' contact is ignored , default 1000, means
#'        for each locus, only consider contacts that within 1000 bins distance away of it
#' @param t1thr Lower bound for t1 for calling TAD, default 0.5 quantile of test statistics
#'        for TADs, 0.9 for subTADs
#' @param fcthr Lower bound quantile for fold change of proportion of 1 between
#'        intra and inter domains, default 0.9
#' @param logt Do log-transformation or not, default TRUE
#' @param dom_order Maximum level of hierarchical structures, default 2 (call TADs and subTADs)
#'         (if all between-within p-value < sthr (after bonferoni correction), stop call TAD or subTAD)
#'   Max_d, Max_dp, Bg_d should be specified in number of bins and change with resloution
#' @return A list includes following elements:
#' \item{tads}{A data frame with columns start, end indicates the start and end coordinates of each domain, respectively}
#' \item{hierTads}{A data frame with columns start, end, dom_order, where dom_order indicates the hierarchical status of a domain, 1 indicates tads, 2 indicates subtads, and so on}
#' \item{params}{A data frame gives the final parameters for calling TADs}
#' @rdname rGMAP
#' @export
rGMAP <- function(hic_mat, resl = 10*10^3, logt = T, dom_order = 2,
                  min_d = 50, Max_d = 200, min_dp = 5, Max_dp = 10,
                  Bg_d = 200, hthr = 0.9,
                  bthr = 400, t1thr = 0.75, fcthr = 0.9){

  if(ncol(hic_mat) == 3){
    names(hic_mat) = c('n1', 'n2', 'counts')


    hic_mat = sparseMatrix(i = hic_mat$n1, j = hic_mat$n2, x=hic_mat$counts,
                           symmetric = T )
  }

  hic_mat = as.matrix(hic_mat)
  hic_mat[is.na(hic_mat)] = 0L
  hic_mat = round(hic_mat)
  options(scipen = 10)


  # remove first k rows and columns if the first k by k submatrix is zero matrix
  k = 0
  for(i in 1:nrow(hic_mat)){
    if(any(hic_mat[i, ] != 0)) {
      k = i - 1
      break
    }
  }
  if(k > 0) {
    hic_mat = hic_mat[-(1:k), -(1:k)]
    #message(paste(">>>> remove first ", k, " bins whose counts are all zero"))
  }
  nn = nrow(hic_mat)

  ## adjust distance effect



  ## do logtransformation
  if(logt){
    hic_mat[which(hic_mat < 0, arr.ind = T)] = 0
    hic_mat = log(hic_mat + 1)
  }


  message(">>>> call TADs...")
  res = call_domain(hic_mat, Max_d, min_d, Max_dp, min_dp, Bg_d, fcthr, hthr,
                    bthr, t1thr)
  tads = res$tads
  if(is.null(tads)) {
    message(">>> No tads: probably because Bg_d is specified too small or bthr too big!")
    message(">>> first try: decrease bthr")
    res = call_domain(hic_mat, Max_d, min_d, Max_dp, min_dp, Bg_d, fcthr, hthr,
                      floor(bthr*0.8), t1thr)
    tads = res$tads
    if(is.null(tads)){
      message(">>> second try: incease Bg_d")
      res = call_domain(hic_mat, Max_d, min_d, Max_dp, min_dp, Bg_d * 2 ,
                        fcthr, hthr,
                        bthr, t1thr)
      tads = res$tads
    }
    if(is.null(tads)){
      message(">>> third try: decrease bthr and double Bg_d")
      res = call_domain(hic_mat, Max_d, min_d, Max_dp, min_dp, Bg_d *2,
                        fcthr, hthr,
                        floor(bthr*0.8), t1thr)
      tads = res$tads
    }

    if(is.null(tads)){
      message(">>> Give up, No tads detected!")
      return(list("tads" = NULL, 'hierTads' = NULL))
    }

  }
  params = data.frame('score' = round(res$score, 1), 'd' = res$d, 'dp' = res$dp, 't1' = res$t1, 't2' = res$t2)

  message(paste(c('score', 'd =', 'dp =', 't1 =', 't2 ='), unlist(params),  collapse = ", "))

  ## call hierarchical structure

  hierTads = list()
  if(dom_order > 1){
    ll = 1
    parenTads = tads
    if(is.null(tads)) return(list("tads" = NULL, 'hierTads' = NULL))


    hierTads[[ll]] = data.frame(tads)

    while(ll < dom_order){
      tempTads = list()
      ll = ll + 1
      for(i in 1:nrow(parenTads)){

        Sbin = parenTads[i, 1]
        Ebin = parenTads[i, 2]
        len = Ebin - Sbin + 1
        if(len <= min(40, floor(500*1000/resl))) next

        message(paste(">>>> call sub-domains for ", i, "th"," domain of order ", ll-1,
                      "...", sep = ""))
        message(paste("start from bin ", parenTads[i, 1], " to bin ", parenTads[i, 2], sep = ""))


        md = min(200*1000/resl, 15)
        Md = min(5 * md, floor(len/3))
        Md = max(md, Md)

        mdp = min_dp
        Mdp = Max_dp

        tmp0 <-  call_domain(hic_mat[Sbin:Ebin, Sbin:Ebin], Md, md, Mdp, mdp,
                             Bg_d = len, fcthr = max(fcthr, 0.925), hthr = max(hthr, 0.99),  len,
                             t1thr = max(0.925, t1thr))

        if(is.null(tmp0$tads)){
          message(paste('>>> no sub-TADs found!'))
          next
        }

        tempTads[[i]] = tmp0$tads + Sbin - 1
        params0 = c(tmp0$d, tmp0$dp, tmp0$t1, tmp0$t2)

        message(paste(c('d =', 'dp =', 't1 =', 't2 ='), params0,  collapse = ", "))

      }

      ## do not go to next step if no subTads are called from current level
      if(length(tempTads) == 0) break

      parenTads = data.frame(do.call("rbind", tempTads), row.names = NULL)

      hierTads[[ll]] = parenTads

    }
    rm(hic_mat)
  }

  message(">>>> Domain calling all done!")

  tads = data.frame((tads+k) * resl - floor(resl/2))

  if(length(hierTads) == 0) {
    hierTads = tads
    hierTads$dom_order = 1
  }else{
    hierTads = lapply(hierTads, function(x) (x + k) * resl - floor(resl/2))

    dlens = sapply(hierTads, nrow)

    hierTads = do.call('rbind', hierTads)

    hierTads$dom_order = rep(1:length(dlens), dlens)
  }


  row.names(hierTads) = NULL
  return(list('tads' = tads, 'hierTads' = hierTads, 'params' = params))

}

rGMAP = cmpfun(rGMAP)





## generate simulated hic_mat and true tads
#' @export
data_simu <- function(stype = 'poisson', nratio = 2.5, mu0 = 20, resl = 1){

  if(stype == 'poisson-dist'){
    ## TADs with gap and Hier
    tbins = 1000  ## total bins

    bounds = c(1, 120, 160, 215, 355, 440, 530, 705, 765, 850, 950, 1000)
    bounds = sort(unique(bounds))
    #mu0 = 5
    ntad = length(bounds) - 1


    n = max(bounds)
    hic_mat  = matrix(0, n, n)
    hic_mat = as.data.table(as.data.frame(hic_mat))

    generateSparse_dat <- function(n1, mu1){
      df = data.table('id1' = rep(1:(n1), (n1):1),
                      'id2' = do.call('c', sapply(1:(n1), function(x) (x):n1)))
      df[, 'dist' := (id2 - id1)]
      df[, dist := ifelse(dist == 0, 0.9, dist)]
      df[, 'N' := rpois(1, mu1  * dist^(-1)), by = list(id1, id2)]

      return(df)
    }

    df = generateSparse_dat(n, mu0)

    for(k in 1:nrow(df)) set(hic_mat, i = as.integer(df$id1[k]),
                             j = as.integer(df$id2[k]), value = df$N[k])
    for(k in 1:nrow(df)) set(hic_mat, i = as.integer(df$id2[k]),
                             j = as.integer(df$id1[k]), value = df$N[k])

    hic_mat = as.matrix(hic_mat)


    modify_dat <- function(hic_mat, bounds, skip_ind = 0, mu = 5){
      ntad = length(bounds) - 1
      sizes = diff(bounds) + 1
      for(i in 1:ntad){
        if(i %in% skip_ind) next
        tn = sizes[i]
        #t_mu = mu * (1 + min(sizes)/tn)
        t_mu = mu
        df = generateSparse_dat(tn, t_mu)
        thic_mat  = matrix(0, tn, tn)
        thic_mat = as.data.table(as.data.frame(thic_mat))
        for(k in 1:nrow(df)) set(thic_mat, i = as.integer(df$id1[k]),
                                 j = as.integer(df$id2[k]), value = df$N[k])
        for(k in 1:nrow(df)) set(thic_mat, i = as.integer(df$id2[k]),
                                 j = as.integer(df$id1[k]), value = df$N[k])
        thic_mat = as.matrix(thic_mat)
        start = bounds[i]
        end = bounds[i+1]
        hic_mat[start:end, start:end] = thic_mat
      }
      return(hic_mat)

    }

    ## add TADs
    hic_mat = modify_dat(hic_mat, bounds, skip_ind = c(5, 7), mu0 * nratio)


    tads_true <- data.frame('start' = bounds[1:ntad], 'end' = bounds[2:(ntad + 1)])
    tads_true <- tads_true[-c(5, 7), ] * resl ## a gap

  }

  if(stype == 'poisson-dist-hier'){
    ## TADs with gap and Hier
    tbins = 1000  ## total bins

    bounds = c(1, 120, 160, 215, 355, 440, 530, 705, 765, 850, 950, 1000)
    bounds = sort(unique(bounds))
    #mu0 = 5
    ntad = length(bounds) - 1


    n = max(bounds)
    hic_mat  = matrix(0, n, n)
    hic_mat = as.data.table(as.data.frame(hic_mat))

    generateSparse_dat <- function(n1, mu1){
      df = data.table('id1' = rep(1:(n1), (n1):1),
                      'id2' = do.call('c', sapply(1:(n1), function(x) (x):n1)))
      df[, 'dist' := (id2 - id1)]
      df[, dist := ifelse(dist == 0, 0.9, dist)]
      df[, 'N' := rpois(1, mu1  * dist^(-1)), by = list(id1, id2)]

      return(df)
    }

    df = generateSparse_dat(n, mu0)

    for(k in 1:nrow(df)) set(hic_mat, i = as.integer(df$id1[k]),
                             j = as.integer(df$id2[k]), value = df$N[k])
    for(k in 1:nrow(df)) set(hic_mat, i = as.integer(df$id2[k]),
                             j = as.integer(df$id1[k]), value = df$N[k])

    hic_mat = as.matrix(hic_mat)


    modify_dat <- function(hic_mat, bounds, skip_ind = 0, mu = 5){
      ntad = length(bounds) - 1
      sizes = diff(bounds) + 1
      for(i in 1:ntad){
        if(i == skip_ind) next
        tn = sizes[i]
        #t_mu = mu * (1 + min(sizes)/tn)
        t_mu = mu
        df = generateSparse_dat(tn, t_mu)
        thic_mat  = matrix(0, tn, tn)
        thic_mat = as.data.table(as.data.frame(thic_mat))
        for(k in 1:nrow(df)) set(thic_mat, i = as.integer(df$id1[k]),
                                 j = as.integer(df$id2[k]), value = df$N[k])
        for(k in 1:nrow(df)) set(thic_mat, i = as.integer(df$id2[k]),
                                 j = as.integer(df$id1[k]), value = df$N[k])
        thic_mat = as.matrix(thic_mat)
        start = bounds[i]
        end = bounds[i+1]
        hic_mat[start:end, start:end] = thic_mat
      }
      return(hic_mat)

    }

    ## add TADs
    hic_mat = modify_dat(hic_mat, bounds, skip_ind = 5, mu0 * nratio)



    # generate subTADs
    start = 530
    end = 705

    sbounds = c(start, start + 50, start + 130, end)

    hic_mat = modify_dat(hic_mat, sbounds, skip_ind = 0, mu = mu0 * nratio * nratio)


    tads_true <- data.frame('start' = bounds[1:ntad],
                            'end' = bounds[2:(ntad + 1)]) * resl
    tads_true = tads_true[-5, ]  ## a gap

    sub_tads <- data.frame('start' = sbounds[1:3],
                           'end' = sbounds[2:4]) * resl

    hierTads = rbind(tads_true, sub_tads)
    hierTads$dom_order = c(rep(1, nrow(tads_true)), rep(2, nrow(sub_tads)))

  }

  if(stype == 'poisson-dist-hier2'){
    ## TADs with gap and Hier
    tbins = 1000  ## total bins

    bounds = c(1, 120, 160, 215, 355, 440, 530, 705, 765, 850, 950, 1000)
    bounds = sort(unique(bounds))
    #mu0 = 5
    ntad = length(bounds) - 1


    n = max(bounds)
    hic_mat  = matrix(0, n, n)
    hic_mat = as.data.table(as.data.frame(hic_mat))

    generateSparse_dat <- function(n1, mu1){
      df = data.table('id1' = rep(1:(n1), (n1):1),
                      'id2' = do.call('c', sapply(1:(n1), function(x) (x):n1)))
      df[, 'dist' := (id2 - id1)]
      df[, dist := ifelse(dist == 0, 0.9, dist)]
      df[, 'N' := rpois(1, mu1  * dist^(-1)), by = list(id1, id2)]

      return(df)
    }

    df = generateSparse_dat(n, mu0)

    for(k in 1:nrow(df)) set(hic_mat, i = as.integer(df$id1[k]),
                             j = as.integer(df$id2[k]), value = df$N[k])
    for(k in 1:nrow(df)) set(hic_mat, i = as.integer(df$id2[k]),
                             j = as.integer(df$id1[k]), value = df$N[k])

    hic_mat = as.matrix(hic_mat)


    modify_dat <- function(hic_mat, bounds, skip_ind = 0, mu = 5){
      ntad = length(bounds) - 1
      sizes = diff(bounds) + 1
      for(i in 1:ntad){
        if(i == skip_ind) next
        tn = sizes[i]
        #t_mu = mu * (1 + min(sizes)/tn)
        t_mu = mu
        df = generateSparse_dat(tn, t_mu)
        thic_mat  = matrix(0, tn, tn)
        thic_mat = as.data.table(as.data.frame(thic_mat))
        for(k in 1:nrow(df)) set(thic_mat, i = as.integer(df$id1[k]),
                                 j = as.integer(df$id2[k]), value = df$N[k])
        for(k in 1:nrow(df)) set(thic_mat, i = as.integer(df$id2[k]),
                                 j = as.integer(df$id1[k]), value = df$N[k])
        thic_mat = as.matrix(thic_mat)
        start = bounds[i]
        end = bounds[i+1]
        hic_mat[start:end, start:end] = thic_mat
      }
      return(hic_mat)

    }

    ## add TADs
    hic_mat = modify_dat(hic_mat, bounds, skip_ind = 5, mu0 * nratio)



    # generate subTADs
    start = 530
    end = 705

    sbounds = c(start, start + 25, end)

    hic_mat = modify_dat(hic_mat, sbounds, skip_ind = 0, mu = mu0 * nratio * nratio)


    tads_true <- data.frame('start' = bounds[1:ntad],
                            'end' = bounds[2:(ntad + 1)]) * resl
    tads_true = tads_true[-5, ]  ## a gap

    sub_tads <- data.frame('start' = sbounds[1:2],
                           'end' = sbounds[2:3]) * resl

    hierTads = rbind(tads_true, sub_tads)
    hierTads$dom_order = c(rep(1, nrow(tads_true)), rep(2, nrow(sub_tads)))

  }


  if(stype == 'nb-dist-hier'){
    ## TADs with gap and Hier
    tbins = 1000  ## total bins

    bounds = c(1, 120, 160, 215, 355, 440, 530, 705, 765, 850, 950, 1000)
    bounds = sort(unique(bounds))
    #mu0 = 5
    ntad = length(bounds) - 1


    n = max(bounds)
    hic_mat  = matrix(0, n, n)
    hic_mat = as.data.table(as.data.frame(hic_mat))

    phi0 = 4
    generateSparse_dat <- function(n1, mu1){
      df = data.table('id1' = rep(1:(n1), (n1):1),
                      'id2' = do.call('c', sapply(1:(n1), function(x) (x):n1)))
      df[, 'dist' := (id2 - id1)]
      df[, dist := ifelse(dist == 0, 0.9, dist)]
      df[, 'N' := rnbinom(1, mu = mu1  * dist^(-1), size = phi0 * mu1  * dist^(-1)), by = list(id1, id2)]

      return(df)
    }

    df = generateSparse_dat(n, mu0)

    for(k in 1:nrow(df)) set(hic_mat, i = as.integer(df$id1[k]),
                             j = as.integer(df$id2[k]), value = df$N[k])
    for(k in 1:nrow(df)) set(hic_mat, i = as.integer(df$id2[k]),
                             j = as.integer(df$id1[k]), value = df$N[k])

    hic_mat = as.matrix(hic_mat)


    modify_dat <- function(hic_mat, bounds, skip_ind = 0, mu = 5){
      ntad = length(bounds) - 1
      sizes = diff(bounds) + 1
      for(i in 1:ntad){
        if(i == skip_ind) next
        tn = sizes[i]
        #t_mu = mu * (1 + min(sizes)/tn)
        t_mu = mu
        df = generateSparse_dat(tn, t_mu)
        thic_mat  = matrix(0, tn, tn)
        thic_mat = as.data.table(as.data.frame(thic_mat))
        for(k in 1:nrow(df)) set(thic_mat, i = as.integer(df$id1[k]),
                                 j = as.integer(df$id2[k]), value = df$N[k])
        for(k in 1:nrow(df)) set(thic_mat, i = as.integer(df$id2[k]),
                                 j = as.integer(df$id1[k]), value = df$N[k])
        thic_mat = as.matrix(thic_mat)
        start = bounds[i]
        end = bounds[i+1]
        hic_mat[start:end, start:end] = thic_mat
      }
      return(hic_mat)

    }

    ## add TADs
    hic_mat = modify_dat(hic_mat, bounds, skip_ind = 5, mu0 * nratio)



    # generate subTADs
    start = 530
    end = 705

    sbounds = c(start, start + 50, start + 130, end)

    hic_mat = modify_dat(hic_mat, sbounds, skip_ind = 0, mu = mu0 * nratio * nratio)


    tads_true <- data.frame('start' = bounds[1:ntad],
                            'end' = bounds[2:(ntad + 1)]) * resl
    tads_true <- tads_true[-5, ]  ## a gap

    sub_tads <- data.frame('start' = sbounds[1:3],
                           'end' = sbounds[2:4]) * resl

    hierTads <- rbind(tads_true, sub_tads)
    hierTads$dom_order = c(rep(1, nrow(tads_true)), rep(2, nrow(sub_tads)))

  }

  if(stype == 'nb-dist'){
    ## TADs with gap and Hier
    tbins = 1000  ## total bins

    bounds = c(1, 120, 160, 215, 355, 440, 530, 705, 765, 850, 950, 1000)
    bounds = sort(unique(bounds))
    #mu0 = 5
    ntad = length(bounds) - 1


    n = max(bounds)
    hic_mat  = matrix(0, n, n)
    hic_mat = as.data.table(as.data.frame(hic_mat))

    phi0 = 4
    generateSparse_dat <- function(n1, mu1){
      df = data.table('id1' = rep(1:(n1), (n1):1),
                      'id2' = do.call('c', sapply(1:(n1), function(x) (x):n1)))
      df[, 'dist' := (id2 - id1)]
      df[, dist := ifelse(dist == 0, 0.9, dist)]
      df[, 'N' := rnbinom(1, mu = mu1  * dist^(-1), size = phi0 * mu1  * dist^(-1)), by = list(id1, id2)]

      return(df)
    }

    df = generateSparse_dat(n, mu0)

    for(k in 1:nrow(df)) set(hic_mat, i = as.integer(df$id1[k]),
                             j = as.integer(df$id2[k]), value = df$N[k])
    for(k in 1:nrow(df)) set(hic_mat, i = as.integer(df$id2[k]),
                             j = as.integer(df$id1[k]), value = df$N[k])

    hic_mat = as.matrix(hic_mat)


    modify_dat <- function(hic_mat, bounds, skip_ind = 0, mu = 5){
      ntad = length(bounds) - 1
      sizes = diff(bounds) + 1
      for(i in 1:ntad){
        if(i %in% skip_ind) next
        tn = sizes[i]
        #t_mu = mu * (1 + min(sizes)/tn)
        t_mu = mu
        df = generateSparse_dat(tn, t_mu)
        thic_mat  = matrix(0, tn, tn)
        thic_mat = as.data.table(as.data.frame(thic_mat))
        for(k in 1:nrow(df)) set(thic_mat, i = as.integer(df$id1[k]),
                                 j = as.integer(df$id2[k]), value = df$N[k])
        for(k in 1:nrow(df)) set(thic_mat, i = as.integer(df$id2[k]),
                                 j = as.integer(df$id1[k]), value = df$N[k])
        thic_mat = as.matrix(thic_mat)
        start = bounds[i]
        end = bounds[i+1]
        hic_mat[start:end, start:end] = thic_mat
      }
      return(hic_mat)

    }

    ## add TADs
    hic_mat = modify_dat(hic_mat, bounds, skip_ind = c(5, 7), mu0 * nratio)


    tads_true <- data.frame('start' = bounds[1:ntad],
                            'end' = bounds[2:(ntad + 1)]) * resl
    tads_true <- tads_true[-c(5, 7), ]  ## a gap


  }

  if(!grepl(stype, pattern = 'hier'))  {
    hierTads = tads_true
    hierTads$dom_order = 1
  }

  row.names(tads_true) = row.names(hierTads) = NULL
  return(list("hic_mat" = hic_mat, "tads_true" = tads_true, 'hierTads' = hierTads))
}
data_simu = cmpfun(data_simu)
