cmp_dist <- function(path_gen = NULL, cmp_dfm = NULL, val = 5, is_impute = TRUE, tEnd = 10) {
  # Generate path from func, compare the distribution with cmp_dfm when it is provided at val
  # returns running time, effective sample size and effective sample size per second
  
  time_total <- system.time(rslt <- path_gen)
  print(time_total)
  time_run <- time_total['elapsed']
  
  if (is_impute == 'Data') {
    rslt <- rslt$xObsMat
  } else if(is_impute) {
    rslt <- timegrid_impute(rslt$dfm, dfm = timegrid(tEnd = tEnd))
    if (is.list(rslt)) {
      rslt <- rslt[[1]]
    }
  }
  
  
  if (is_impute == 'Data') {
    es <- effectiveSize(t(rslt))
  } else if(is_impute) {
    es <- effectiveSize(rslt[rslt$time == val, ]$omega)
  }
  ess <- es / time_run
  
  if (length(nrow(cmp_dfm)) != 0) {
    ks_rslt <- ks.test(rslt[rslt$time == val, ]$omega, 
                       cmp_dfm[cmp_dfm$time == val, ]$omega)
    
    ks_dist <- ks_rslt$statistic
    ks_pval <- ks_rslt$p.value
    
    return(c(time_run, es, ess, ks_dist, ks_pval))
  }
  
  return(c(time_run, es, ess))
}
