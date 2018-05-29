##---------------------------------------------------------------------##
## --This program is written by Wenbao Yu for implementing GMAP -------##
## --This file contains all the source codes
## -- filter peak if it's not a local peak of insulation in orignal; update tune_score
## -- function without using prefixed domain as background
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


  return(list('up' = round(stat_up, 2), 'down' = round(stat_down,2),
              'wb' = round(stat_wb, 2),
              'pr_fc_wb' = round(pr_fc_wb, 2)))

}
cal_stat = cmpfun(cal_stat)



## find local peaks for a given vector
localMaxima <- function(x,  stats1, thr = 0.75,  dp = 10) {
  # smmothing using sliding window before search local peaks
  # filter by fc
  x = round(runmean(x, 5), 2)

  stats1 = round(runmean(stats1, 5), 2)

  y = extrema(x)$maxindex[, 1]

  y_fc = extrema(stats1)$maxindex[, 1]

  y_mdist_fc = sapply(y, function(t) min(abs(t - y_fc)))

  y = y[y_mdist_fc <= 5]
  xy = x[y]

  # filter by t1
  local_t1_quantile <- function(ly){
    ll = max(ly - 500, 1)
    rr = min(ly + 500, length(x))
    return(quantile(x[ll : rr], thr))
  }

  local_t1thr = quantile(x, thr)

  if(length(x) >= 1000) local_t1thr = sapply(y, local_t1_quantile)


  y = y[xy >= local_t1thr]  ## screening by t1


  if(length(y) == 0) return(NULL)

  # remove peaks that are very close
  if(!is.null(y)) y = combPeaks(y, x, dp)


  y = y[y >= 6]
  y = y[y <= (length(x) - 6)]

  return(y)
}
localMaxima = cmpfun(localMaxima)



## define tad by gmap given peaks and statistics
def_tad_gmap <- function(wb_peaks, stat_up, stat_down, nn, thr = 2){
  s1 = stat_up[wb_peaks]
  s2 = stat_down[wb_peaks]


  len = length(wb_peaks)
  signs = rep(0, len)

  ind1 = which(s1 > thr)
  ind2 = which(s2 > thr)
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
  bounds = sort(bounds)
  len = length(bounds)
  return(data.frame('start' = bounds[-len], 'end' = bounds[-1]))
}
full_seg = cmpfun(full_seg)

## calculate insulation score
cal_insul = function(x, pp, wd = 50){
  nn = nrow(pp)
  start_bin = max(1, (x - wd))
  end_bin = min(nn, (x + wd))
  intra.score = sum(c(pp[start_bin:x, start_bin:x], pp[x:end_bin, x:end_bin]))
  inter.score = 2*sum(pp[start_bin:x, x:end_bin])
  insul.score = log2(intra.score + 1) - log2(inter.score + 1)
  return(insul.score)
}


## without using predifined area: using intra-tad over inter-tad instead
tune_insulScore <- function(pp, tads, wd = 50){
  nn = nrow(pp)

  bds = sort(unique(unlist(tads)))
  bds = setdiff(bds, c(1, nn))
  insul.score = lapply(bds, cal_insul, pp, wd)
  insul.score = do.call('c', insul.score)

  if(all(is.na(insul.score))) return(0)

  if(F){
    score.high = insul.score[insul.score > 2]
    score.midd = insul.score[insul.score > 1 & insul.score <= 2]
    score.low = insul.score[insul.score < 0]

    score.summ = 0
    if(length(score.high) > 0) score.summ = score.summ + 2 * mean(score.high) * log2(length(score.high))
    if(length(score.low) > 0) score.summ = score.summ + 2 * mean(score.low) * log2(length(score.low))
    if(length(score.midd) > 0) score.summ = score.summ + mean(score.midd) * log2(length(score.midd))

  }

  score.summ = mean(insul.score) * log(length(insul.score[insul.score > 0]))
  return(score.summ)
}
tune_insulScore = cmpfun(tune_insulScore)

## including tune t1 and t2 -- given dp and d
tune_T1T2 <- function(pp, stats, stats1, d, dp, t1thr = 0.5){

  stat_up = stats$up
  stat_down = stats$down
  stat_wb = stats$wb
  pr_fc_wb = stats$pr_fc_wb
  nn = nrow(pp)

  # no local peaks
  if(all(diff(stat_wb) <=0) || all(diff(stat_wb) >= 0)) return(list('score' = 0, 'tads' = NULL))

  # 75% quantile is too close to 99% quantile
  if(max(stat_wb)- quantile(stat_wb, 0.75) < 1) return(list('score' = 0, 'tads' = NULL))

  sthr = 0.05
  sthr0 = round(qnorm(sthr/nn, lower.tail = F), 1)

  # give candidates for t1
  stat_wb = runmean(stat_wb, 5)

  cand_t1 = seq(0.99, t1thr, length = 50)  # use quantile

  # define t2 candidates globally
  tmp = c(stat_up, stat_down)
  cand_t2 <- quantile(tmp[tmp != 0], seq(0.975, 1.0, length = 2))
  tmp = qnorm(sthr, lower.tail = F)
  cand_t2 = unique(pmax(tmp, round(cand_t2, 2)))


  len1 = length(cand_t1)
  len2 = length(cand_t2)
  score = matrix(0, len1, len2)

  for(i in 1:len1){
    wb_peaks = localMaxima(stat_wb, stats1, thr = cand_t1[i], dp = dp)

    if(is.null(wb_peaks)) {
      score[i, ] = 0
      next
    }

    if(length(wb_peaks) == 0) {
      score[i, ] = 0
      next
    }

    for(j in 1:len2){

      tads <- def_tad_gmap(wb_peaks, stat_up, stat_down, nn, thr = cand_t2[j])

      if(nrow(tads) == 0) {
        score[i, j] = 0
      }else{
        score[i, j] <- tune_insulScore(pp, tads)
      }
    }
  }

  if(max(score, na.rm = TRUE) <= 0) return(list('score' = 0, 'tads' = NULL))

  ind = which(score == max(score, na.rm = TRUE), arr.ind=TRUE)[1, ]
  wb_peaks = localMaxima(stat_wb, stats1, thr = cand_t1[ind[1]],  dp = dp)

  if(is.null(wb_peaks)) return(list('score' = 0, 'tads' = NULL))

  tads <- def_tad_gmap(wb_peaks, stat_up, stat_down, nn,
                       thr = cand_t2[ind[2]])

  return(list('tads' = tads, 'score' = max(score), 't1' = round(cand_t1[ind[1]], 3),
              't2' = round(cand_t2[ind[2]], 3) , 'stats' = stats,
              'local_peaks' = wb_peaks))
}
tune_T1T2 = cmpfun(tune_T1T2)




## tune all parameters -- give candidates of d and dp
tune_allPara <- function(pp, pp1, candD, candDp,  t1thr = 0.5){
  lend = length(candD)
  lendp = length(candDp)
  score = matrix(0, lend, lendp)
  nn = nrow(pp)
  sthr = 0.05

  for(i in 1:lend){
    di = candD[i]
    stats = cal_stat(pp, d = di)
    stats1 = lapply(1:nn, cal_insul, pp1, 25)
    stats1 = do.call('c', stats1)

    if(all(stats$wb <= qnorm(sthr/nn, lower.tail = F))) {
      score[i, ] = 0
      next
    }
    for(j in 1:lendp){
      score[i, j] = tune_T1T2(pp, stats, stats1, d = candD[i], dp = candDp[j],
                               t1thr = t1thr)$score
    }
  }

  if(max(score, na.rm = TRUE) <= 0) return(NULL)

  ind = which(score == max(score, na.rm = TRUE), arr.ind = T)[1, ]
  d = candD[ind[1]]
  dp = candDp[ind[2]]

  stats = cal_stat(pp, d = d)
  res = tune_T1T2(pp, stats,  stats1, d = d, dp = dp, t1thr = t1thr)
  res$d = d
  res$dp = dp
  return(res)

}
tune_allPara = cmpfun(tune_allPara)



## remove redundant domains
rm_fpdomain <- function(subTads, pp){

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

  ids = which(s_density <= e_density)
  if(length(ids) > 0) subTads = subTads[-ids, ]

  return(subTads)
}
rm_fpdomain = cmpfun(rm_fpdomain)


## remove redundant domains
## dw represents dw up- and down- range of the current domain
rm_fpdomain_local <- function(subTads, pp, dw = 5){

  subTads = as.matrix(subTads)
  bds = sort(unique(as.vector(subTads)))

  ## filter domains with very small contacts
  subTads = as.matrix(full_seg(bds))

  s_density = round(apply(subTads, 1, function(x) mean(pp[x[1]:x[2], x[1]:x[2]])), 2)
  tsize = round(apply(subTads, 1, function(x) x[2] - x[1]))
  midp = apply(subTads, 1, function(x) floor(x[2]/2 + x[1]/2))
  nn = nrow(pp)

  calDensity_expect_local <- function(j, dw0 = dw){
    tsize0 = tsize[j]
    midp0 = midp[j]
    ids = lapply(max(1, midp0 - dw0 * tsize0) : min(nn, midp0 + dw0 * tsize0), function(x) return(cbind(x, c(max(1, x - tsize0):min(nn, x + tsize0)))))
    ids = do.call('rbind', ids)
    return(round(mean(pp[ids]), 2))
  }


  e_density = sapply(1:length(tsize), calDensity_expect_local)

  ids = which(s_density <= e_density)
  if(length(ids) > 0) subTads = subTads[-ids, ]

  return(subTads)
}
rm_fpdomain_local = cmpfun(rm_fpdomain_local)


## transfrom hic normalized count to binary
gmixmodel <- function(sub_mat, hthr = 0.9, bthr = 200){
  nn = nrow(sub_mat)
  pp = matrix(0L, nn, nn)

  # if the tad is too small (less than 10 bins), do not call subtad
  if(nn <= 10) return(pp)

    ## if the hic-matrix is too large, we assume two loci with distance larger than
    ## say 1000 bins have no interactions. Data from this part will not be used for
    ## constructing the mixture models
    if(all(class(sub_mat) != "data.frame")) {
      mat = data.frame(sub_mat)
    }
    mat = as.list(mat)  ## for faster subsetting

    colIDs =  sapply(1:(nn - 1), function(x) return(c((x + 1) : min(nn, x + bthr))))
    len = sapply(colIDs, length)

    getData <- function(x){
      return(mat[[x]][colIDs[[x]]])
    }
    temp = unlist(sapply(1:(nn - 1), getData))
    rm(mat)



    rowIDs = rep(1:(nn - 1), len)
    colIDs = unlist(colIDs)

    partialIDs = cbind(rowIDs, colIDs)
    rm(rowIDs, colIDs)

    # deal with outliers
    class1 = rep(0, length(temp))

    # deal with outliers

    out.thr.upper = quantile(temp, 0.975)
    out.thr.lower = max(quantile(temp, 0.025), 0)

    id.lower = which(temp <= out.thr.lower)
    id.upper = which(temp >= out.thr.upper)
    ids = NULL

    if(length(id.upper) > 0){
      class1[id.upper] = 1
      ids = c(id.lower, id.upper)
      temp = temp[-ids]
      class0 = class1[-ids]
    }else{
      if(length(id.lower) > 0){
        ids = id.lower
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

    if(length(ids) > 0) {
      class1[-ids] = class0
      pp[partialIDs] = class1
    }else{
      pp[partialIDs] = class0
    }

    pp = pp + t(pp)
    rm(class0, partialIDs)



  diag(pp) <- 0

  return(pp)
}
gmixmodel = cmpfun(gmixmodel)


## use default classification; using hthr quantile to cluster
## only when gaussian mixture model is not converged
gmixmodel0 <- function(sub_mat, hthr = 0.9, bthr = 200){
  nn = nrow(sub_mat)
  pp = matrix(0L, nn, nn)

  # if the tad is too small (less than 10 bins), do not call subtad
  if(nn <= 10) return(pp)

  # do mixture model clustering

    ## if the hic-matrix is too large, we assume two loci with distance larger than
    ## say 1000 bins have no interactions. Data from this part will not be used for
    ## constructing the mixture models
    if(all(class(sub_mat) != "data.frame")) {
      mat = data.frame(sub_mat)
    }
    mat = as.list(mat)  ## for faster subsetting

    colIDs =  sapply(1:(nn - 1), function(x) return(c((x + 1) : min(nn, x + bthr))))
    len = sapply(colIDs, length)

    getData <- function(x){
      return(mat[[x]][colIDs[[x]]])
    }
    temp = unlist(sapply(1:(nn - 1), getData))
    rm(mat)



    rowIDs = rep(1:(nn - 1), len)
    colIDs = unlist(colIDs)

    partialIDs = cbind(rowIDs, colIDs)
    rm(rowIDs, colIDs)

    # deal with outliers
    class1 = rep(0, length(temp))

    # deal with outliers

    out.thr.upper = quantile(temp, 0.975)
    out.thr.lower = max(quantile(temp, 0.025), 0)

    id.lower = which(temp <= out.thr.lower)
    id.upper = which(temp >= out.thr.upper)
    ids = NULL

    if(length(id.upper) > 0){
      class1[id.upper] = 1
      ids = c(id.lower, id.upper)
      temp = temp[-ids]
      class0 = class1[-ids]
    }else{
      if(length(id.lower) > 0){
        ids = id.lower
        temp = temp[-ids]
        class0 = class1[-ids]
      }else{
        class0 = class1
      }
    }


    if(length(unique(temp)) >= 2){
      model2 <- Mclust(temp, G = 2, modelNames = 'E')

      if(!is.null(model2)){
        class0 = model2$classification - 1
        if(mean(temp[class0 == 1]) <= mean(temp)) class0 = 1 - class0
      }

      if(is.null(model2)){ ## not covergent
        message('Trying using models with different variances')
        model2 <- Mclust(temp, G = 2, modelNames = 'V')
        if(!is.null(model2)) {
          class0 = model2$classification - 1
          if(mean(temp[class0 == 1]) <= mean(temp)) class0 = 1 - class0
        }else{
          message('No models can be fitted, using hthr quantile as cutoff')
          class0 = ifelse(temp > quantile(temp, hthr), 1, 0)
        }
      }

      # just in case:
      if(length(unique(class0)) == 1 )  class0 = ifelse(temp > quantile(temp, hthr), 1, 0)
      rm(temp, model2)
    }

    if(length(ids) > 0) {
      class1[-ids] = class0
      pp[partialIDs] = class1
    }else{
      pp[partialIDs] = class0
    }

    pp = pp + t(pp)
    rm(class0, partialIDs)

  diag(pp) <- 0

  return(pp)
}
gmixmodel0 = cmpfun(gmixmodel0)


call_domain <- function(sub_mat, Max_d, min_d, Max_dp, min_dp,  hthr = 0.5,
                        bthr = 200, t1thr = 0.5){

  nn = nrow(sub_mat)
  ##
  pp = gmixmodel(sub_mat, hthr, bthr)

  #write.table(pp, file = 'tmp/tmp2.txt', row.names = F, col.names = F, quote = F, sep = '\t')

  if(sum(pp) == 0) return(list('tads' = NULL))


  # define candidates for d and dp
  #candD = seq(min_d, Max_d, by = min_d)
  candD = seq(min_d, Max_d, by = min_d)

  candDp = seq(min_dp, Max_dp, by = min_dp)

  # define the background
  if(bthr > nn){

    rnum_bk <- sum(pp)
    tnum_bk <- length(pp)

    rnum_bk0 = sum(sub_mat)

  }else{
    idBg = lapply(1:nn, function(x) return(cbind(x, c(max(1, x - bthr):min(nn, x + bthr)))))
    idBg = do.call('rbind', idBg)
    temp = pp[idBg]
    rnum_bk <- sum(temp)
    tnum_bk <- length(temp)
  }

  if(rnum_bk >= 0.95 * tnum_bk || rnum_bk < 1/sqrt(tnum_bk) * tnum_bk) return(list("tads" = NULL))

  pp1 = pp
  pp1[idBg] = sub_mat[idBg]
  diag(pp1) = 0

  res = tune_allPara(pp, pp1, candD, candDp, t1thr)
  subTads = res$tads



  if(!is.null(subTads)) res$tads = rm_fpdomain_local(subTads, pp)

  return(res)
}
call_domain = cmpfun(call_domain)





#' Detect hierarchical choromotin domains by GMAP
#' @param  hic_mat n by n matrix, n is the number of bins. Or hic_mat could by 3 columns
#' matrix or data.frame with columns: bin1, bin2, counts, in which bin1 and bin2, from 1 to m, are the bin
#' index of a hic contact
#' @param  resl The resolution (bin size), default 10kb
#' @param  min_d The minimum d (d: window size), default 25
#' @param  Max_d The maximum d (d: window size), default 100
#' @param min_dp The minmum dp (dp: lower bound of tad size), defalt 5
#' @param Max_dp The maximum dp (dp: lower bound of tad size), defalt 10
#' @param hthr The lower bound cutoff for posterior probability, default 0.99
#' @param bthr Threshold, above which a locus' contact is ignored , default 1000, means
#'        for each locus, only consider contacts that within 1000 bins distance away of it
#' @param t1thr Lower bound for t1 for calling TAD, default 0.25 quantile of test statistics
#'        for TADs, 0.9 for subTADs
#' @param logt Do log-transformation or not, default TRUE
#' @param dom_order Maximum level of hierarchical structures, default 2 (call TADs and subTADs)
#'         (if all between-within p-value < sthr (after bonferoni correction), stop call TAD or subTAD)
#'   Max_d, and Max_dp should be specified in number of bins
#' @return A list includes following elements:
#' \item{tads}{A data frame with columns start, end indicates the start and end coordinates of each domain, respectively}
#' \item{hierTads}{A data frame with columns start, end, dom_order, where dom_order indicates the hierarchical status of a domain, 1 indicates tads, 2 indicates subtads, and so on}
#' \item{params}{A data frame gives the final parameters for calling TADs}
#' @rdname rGMAP
#' @export
rGMAP <- function(hic_mat, resl = 10*10^3, logt = T, dom_order = 2,
                  min_d = 25, Max_d = 100, min_dp = 5, Max_dp = 10,
                  hthr = 0.99, bthr = min(400, 400*10^3/resl), t1thr = 0.25){

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

  ## do logtransformation
  if(logt){
    hic_mat[which(hic_mat < 0, arr.ind = T)] = 0
    hic_mat = log2(hic_mat + 1)
  }

  message(">>>> call TADs...")
  res = call_domain(hic_mat, Max_d, min_d, Max_dp, min_dp, hthr, bthr, t1thr)
  tads = res$tads
  if(is.null(tads)) {
    message(">>> No tads: probably because bthr, t1thr too big!")

    message(">>> first try: decrease bthr")
    res = call_domain(hic_mat, Max_d, min_d, Max_dp, min_dp, hthr,
                      floor(bthr*0.5), t1thr)
    tads = res$tads

    if(is.null(tads)){
      message(">>> second try: decrease t1thr")
      res = call_domain(hic_mat, Max_d, min_d, Max_dp, min_dp, hthr,
                        bthr, t1thr = t1thr/2)
      tads = res$tads
    }

    if(is.null(tads)){
      message(">>> third try: decrease bthr and t1thr")
      res = call_domain(hic_mat, Max_d, min_d, Max_dp, min_dp,  hthr,
                        floor(bthr*0.5), t1thr = t1thr/2)
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


        md = min(200*1000/resl, min_d)
        Md = min(Max_d, floor(len/3))
        Md = max(md, Md)

        mdp = min_dp
        Mdp = Max_dp

        tmp0 <-  call_domain(hic_mat[Sbin:Ebin, Sbin:Ebin], Md, md, Mdp, mdp,
                             hthr = max(hthr, 0.9), floor(len/2), t1thr = max(t1thr, 0.9))

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



## transfer tad to be in triangle shape
triangle_tad <- function(tads, resl = 1){
  pos <- within(tads, {
    y <- end/resl
    x <- start/resl
    rm(start, end)
  })

  # to plot a tad in a triangle shape
  triangle_tran <- function(x){
    yy = matrix(0, 3, 2)
    yy[1, ] = c(x[1], x[1])
    yy[2, ] = c(x[1], x[2])
    yy[3, ] = c(x[2], x[2])
    return(yy)
  }

  xy = NULL
  pos = as.matrix(pos)
  for(i in 1:nrow(pos)){
    xy = rbind(xy, triangle_tran(pos[i, ]))
  }
  return(xy)
}


## transferd triangle TADs for better visualization (horizontally arranged)
transf_tad_horz <- function(triangle_tad){
  xy_tad = as.matrix(triangle_tad)
  transf_mat = matrix(c(1, -1, 1, 1), 2, 2)/2
  #transf_mat = matrix(c(1, -1, 1, 1), 2, 2)/sqrt(2)
  xy_new =  transf_mat %*% t(xy_tad)
  return(t(xy_new))
}




#' visualize hierarchical domains
#' @param  hic_mat hic_mat with 3 columns
#' matrix or data.frame with columns: bin1, bin2, counts, in which bin1 and bin2, from 1 to m, are the bin
#' @param  hiertads_gmap the hierarchical domains called by gmap
#' @param start_bin the start bin of the genome
#' @param end_bin the end bin of the genome
#' @param cthr the count threshold for color, default 20
#' @param kb_resl reslution of Hi-C data in kb
#' @rdname plotdom
#' @export
plotdom <- function(hic_dat, hiertads_gmap, start_bin, end_bin, cthr = 20, kb_resl = 10){

  resl = kb_resl * 1000

  names(hiertads_gmap) = c('start', 'end', 'dom_order')
  names(hic_dat) = c('n1', 'n2', 'count')


  tads_gmap = subset(hiertads_gmap, dom_order == 1, select = -dom_order)
  tadsL2 = subset(hiertads_gmap, dom_order == 2, select = -dom_order)
  tadsL3 = subset(hiertads_gmap, dom_order == 3, select = -dom_order)

  ## plot tads_gmap (select region)
  xy_gmap = triangle_tad(tads_gmap, resl = resl)
  if(nrow(tadsL2) > 0 )xy_tadsL2 = triangle_tad(tadsL2, resl = resl)
  if(nrow(tadsL3) > 0 ) xy_tadsL3 = triangle_tad(tadsL3, resl = resl)


  ## plot it horizontally
  hxy_gmap = transf_tad_horz(xy_gmap)
  if(nrow(tadsL2) > 0 ) hxy_tadsL2 = transf_tad_horz(xy_tadsL2)
  if(nrow(tadsL3) > 0 )hxy_tadsL3 = transf_tad_horz(xy_tadsL3)


  mm = start_bin
  MM = end_bin
  dat = data.frame(hic_dat)
  pdat = subset(dat, n1<=MM & n1 >= mm & n2 <= MM & n2 >= mm)
  pdat$count = ifelse(pdat$count >cthr, cthr, pdat$count)

  orgPlot <- ggplot(data = pdat, aes(n1, n2)) + geom_point(aes(colour = count)) +
    scale_colour_gradient(high='red', low = 'white') + xlim(min(pdat$n1), max(pdat$n1))

  #orgPlot


  ## plot it in horizontal triangle
  hdat = pdat[pdat$n1 <= pdat$n2, ]
  hdat[, 1:2] = transf_tad_horz(as.matrix(hdat[, 1:2]))

  hdat = hdat[hdat$n2 <= 100, ]
  #hdat$count = ifelse(hdat$count > 30, 30, hdat$count)


  cols = colorRampPalette(c("white", "red"))(2)

  #library(RColorBrewer)
  #cols <- brewer.pal(10, 'Reds')
  ss = 10^6/resl       ## plot in scale of mb
  hdat$n1 = hdat$n1/ss
  #dd = sample(1:nrow(hdat), 20000)
  #hdat = hdat[dd, ]
  orgPlot <- ggplot(data = hdat, aes(n1, n2)) + geom_point(aes(colour = count)) +
    scale_colour_gradient(high = 'red', low = 'white') + xlim(min(hdat$n1), max(hdat$n1)) +
    xlab('Mb') + ylab("") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                                  legend.position = 'none') + theme(
                                    plot.background = element_blank()
                                    ,panel.grid.major = element_blank()
                                    ,panel.grid.minor = element_blank()
                                    ,panel.border = element_blank())

  #orgPlot



  orgWTad_horz <- function(xy, m = mm, M = MM){
    ind1 = which(xy[, 1]>= m & xy[, 1]<= M)
    if(length(ind1) <= 1) return(NULL)
    xy = xy[ind1, ]
    df = data.frame('xs' = xy[1:(nrow(xy)-1), 1]/ss,
                    'ys' = xy[1:(nrow(xy)-1), 2],
                    "xe" = xy[2:nrow(xy), 1]/ss,
                    "ye" = xy[2:nrow(xy), 2])
    return(df)
  }

  p1 = p2 = p3 = orgPlot

  if(!is.null(orgWTad_horz(hxy_gmap))) {
    tdat = orgWTad_horz(hxy_gmap)
    tdat = tdat[order(tdat$xs), ]
    tdat = tdat[tdat$xs <= tdat$xe,]
    p1 <- p1 + geom_segment(aes(x = xs, y = ys, xend = xe, yend = ye), size = 0.7,
                            data = tdat, color = 'black')
    p2 = p1
  }

  if(nrow(tadsL2) > 0 ){
    if(!is.null(orgWTad_horz(hxy_tadsL2))){
      p2 = p2 + geom_segment(aes(x = xs, y = ys, xend = xe, yend = ye), size = 0.7,
                             data = orgWTad_horz(hxy_tadsL2), color = 'black')
      p3 = p2
    }
  }


  if(nrow(tadsL3) > 0 ){
    if(!is.null(orgWTad_horz(hxy_tadsL3))) {
      p3 = p2 + geom_segment(aes(x = xs, y = ys, xend = xe, yend = ye), size = 0.7,
                             data = orgWTad_horz(hxy_tadsL3), color = 'black')
    }
  }


  return(list('p1' = p1, 'p2' = p2, 'p3' = p3))
}
plotdom = cmpfun(plotdom)

