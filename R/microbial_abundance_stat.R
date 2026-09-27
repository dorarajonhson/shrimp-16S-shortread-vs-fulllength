std.error.of.diff <- function(x1, x2) {
  n1 = length(x1)
  n2 = length(x2)
  m1 = mean(x1)
  m2 = mean(x2)
  sd1 = sd(x1)
  sd2 = sd(x2)
  sem1 = sd(x1)/sqrt(n1)
  sem2 = sd(x2)/sqrt(n2)
  ssq1 = 0; for(x in x1) {ssq1 = ssq1 + x*x}
  ssq2 = 0; for(x in x2) {ssq2 = ssq2 + x*x}
  sed = sqrt(sem1*sem1 + sem2*sem2)
  return(sed)
}

assess.microbe.stat <- function(dat, type) {
  gr1 = dat[dat$pond == 'A' & dat$type == type,]
  gr2 = dat[dat$pond == 'B' & dat$type == type,]
  mean1 = mean(gr1$abundance, na.rm = T)
  mean2 = mean(gr2$abundance, na.rm = T)
  sed = std.error.of.diff(gr1$abundance, gr2$abundance)
  pval = wilcox.test(gr1$abundance, gr2$abundance, exact = F)$p.value
  result = c(mean_1 = mean1, mean_2 = mean2, sed = sed, pvalue = pval)
  return(result)
}