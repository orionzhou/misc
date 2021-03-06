source('functions.R')

dirw = glue('{dird}/16_ase')
#{{{ functions
require(stats4)
require(multidplyr)
require(bbmle)
calc_bic <- function(i, p1, p2, h1, h2, disp, method='L-BFGS-B') {
    #{{{
    cat(i, "\n", sep = '')
    #{{{ LLs
LL1 <- function(mu1, mu3) {
    #{{{
    l1 = dnbinom(p1, size=size, mu=mu1, log=T)
    l2 = dnbinom(p2, size=size, mu=mu1, log=T)
    l3 = dnbinom(h1, size=size, mu=mu3, log=T)
    l4 = dnbinom(h2, size=size, mu=mu3, log=T)
    l = -sum(l1, l2, l3, l4)
    ifelse(is.nan(l), 100, l)
    #}}}
}
LL2 <- function(mu1, mu3, mu4) {
    #{{{
    l1 = dnbinom(p1, size=size, mu=mu1, log=T)
    l2 = dnbinom(p2, size=size, mu=mu1, log=T)
    l3 = dnbinom(h1, size=size, mu=mu3, log=T)
    l4 = dnbinom(h2, size=size, mu=mu4, log=T)
    l = -sum(l1, l2, l3, l4)
    ifelse(is.nan(l), 100, l)
    #}}}
}
LL3 <- function(mu1, mu2, mu3) {
    #{{{
    mu4 = (mu2 / mu1) * mu3
    l1 = dnbinom(p1, size=size, mu=mu1, log=T)
    l2 = dnbinom(p2, size=size, mu=mu2, log=T)
    l3 = dnbinom(h1, size=size, mu=mu3, log=T)
    l4 = dnbinom(h2, size=size, mu=mu4, log=T)
    l = -sum(l1, l2, l3, l4)
    ifelse(is.nan(l), 100, l)
    #}}}
}
LL4 <- function(mu1, mu2, mu3) {
    #{{{
    l1 = dnbinom(p1, size=size, mu=mu1, log=T)
    l2 = dnbinom(p2, size=size, mu=mu2, log=T)
    l3 = dnbinom(h1, size=size, mu=mu3, log=T)
    l4 = dnbinom(h2, size=size, mu=mu3, log=T)
    l = -sum(l1, l2, l3, l4)
    ifelse(is.nan(l), 100, l)
    #}}}
}
LL5 <- function(mu1, mu2, mu3, mu4) {
    #{{{
    l1 = dnbinom(p1, size=size, mu=mu1, log=T)
    l2 = dnbinom(p2, size=size, mu=mu2, log=T)
    l3 = dnbinom(h1, size=size, mu=mu3, log=T)
    l4 = dnbinom(h2, size=size, mu=mu4, log=T)
    l = -sum(l1, l2, l3, l4)
    ifelse(is.nan(l), 100, l)
    #}}}
}
    #}}}
    size = 1 / disp
    n_obs = length(p1)
    m1 = round(mean(p1)); m2 = round(mean(p2)); m3 = round(mean(h1)); m4 = round(mean(h2))
    p1 = round(p1); p2 = round(p2); h1 = round(h1); h2 = round(h2)
    min_mu = 1e-2; max_mu = 1e8
    fit1 = mle2(LL1, start = list(mu1=(m1+m2)/2, mu3=(m3+m4)/2),
          lower = rep(min_mu, 2), upper=rep(max_mu, 2),
          method = method)#, nobs = n_obs)
    fit2 = mle2(LL2, start = list(mu1=(m1+m2)/2, mu3=m3, mu4=m4),
          lower = rep(min_mu, 3), upper=rep(max_mu, 3),
          method = method)#, nobs = n_obs)#, control = list(trace=3, maxit=1000))
    fit3 = mle2(LL3, start = list(mu1=m1, mu2=m2, mu3=(m3+m4)/2),
          lower = rep(min_mu, 3), upper=rep(max_mu, 3),
          method = method)#, nobs = n_obs)
    fit4 = mle2(LL4, start = list(mu1=m1, mu2=m2, mu3=(m3+m4)/2),
          lower = rep(min_mu, 3), upper=rep(max_mu, 3),
          method = method)#, nobs = n_obs)
    fit5 = mle2(LL5, start = list(mu1=m1, mu2=m2, mu3=m3, mu4=m4),
          lower = rep(min_mu, 4), upper=rep(max_mu, 4),
          method = method)#, nobs = n_obs)
    #coef(fitc)
    bic = AICtab(fit1, fit2, fit3, fit4, fit5, k = log(n_obs), sort=F)
    tb = as_tibble(bic) %>%
        mutate(reg = c('conserved','unexpected','cis','trans','cis+trans')) %>%
        arrange(dAIC)
    tb$reg[1]
    #}}}
}
calc_bic_2 <- function(i, p1,p2,h1,h2,bp1,bp2,bh1,bh2,disp, method='L-BFGS-B') {
    #{{{
    cat(i, "\n", sep = '')
    #{{{ LLs
LL1 <- function(mu1,mu2,mu3,mu4,mu5,mu7) {
    #{{{
    mu6 = mu5 - mu1 + mu2
    mu8 = mu7 - mu3 + mu4
    l1 = dnbinom(p1, size=size, mu=mu1, log=T)
    l2 = dnbinom(p2, size=size, mu=mu2, log=T)
    l3 = dnbinom(h1, size=size, mu=mu3, log=T)
    l4 = dnbinom(h2, size=size, mu=mu4, log=T)
    l5 = dnbinom(bp1, size=size, mu=mu5, log=T)
    l6 = dnbinom(bp2, size=size, mu=mu6, log=T)
    l7 = dnbinom(bh1, size=size, mu=mu7, log=T)
    l8 = dnbinom(bh2, size=size, mu=mu8, log=T)
    l = -sum(l1, l2, l3, l4, l5, l6, l7, l8)
    ifelse(is.nan(l), 100, l)
    #}}}
}
LL2 <- function(mu1,mu2,mu3,mu4,mu5,mu7,mu8) {
    #{{{
    mu6 = mu5 - mu1 + mu2
    l1 = dnbinom(p1, size=size, mu=mu1, log=T)
    l2 = dnbinom(p2, size=size, mu=mu2, log=T)
    l3 = dnbinom(h1, size=size, mu=mu3, log=T)
    l4 = dnbinom(h2, size=size, mu=mu4, log=T)
    l5 = dnbinom(bp1, size=size, mu=mu5, log=T)
    l6 = dnbinom(bp2, size=size, mu=mu6, log=T)
    l7 = dnbinom(bh1, size=size, mu=mu7, log=T)
    l8 = dnbinom(bh2, size=size, mu=mu8, log=T)
    l = -sum(l1, l2, l3, l4, l5, l6, l7, l8)
    ifelse(is.nan(l), 100, l)
    #}}}
}
LL3 <- function(mu1,mu2,mu3,mu4,mu5,mu6,mu7) {
    #{{{
    mu8 = (mu6-mu2) / (mu5-mu1) * (mu7-mu3) + mu4
    l1 = dnbinom(p1, size=size, mu=mu1, log=T)
    l2 = dnbinom(p2, size=size, mu=mu2, log=T)
    l3 = dnbinom(h1, size=size, mu=mu3, log=T)
    l4 = dnbinom(h2, size=size, mu=mu4, log=T)
    l5 = dnbinom(bp1, size=size, mu=mu5, log=T)
    l6 = dnbinom(bp2, size=size, mu=mu6, log=T)
    l7 = dnbinom(bh1, size=size, mu=mu7, log=T)
    l8 = dnbinom(bh2, size=size, mu=mu8, log=T)
    l = -sum(l1, l2, l3, l4, l5, l6, l7, l8)
    ifelse(is.nan(l), 100, l)
    #}}}
}
LL4 <- function(mu1,mu2,mu3,mu4,mu5,mu6,mu7) {
    #{{{
    mu8 = mu7 - mu3 + mu4
    l1 = dnbinom(p1, size=size, mu=mu1, log=T)
    l2 = dnbinom(p2, size=size, mu=mu2, log=T)
    l3 = dnbinom(h1, size=size, mu=mu3, log=T)
    l4 = dnbinom(h2, size=size, mu=mu4, log=T)
    l5 = dnbinom(bp1, size=size, mu=mu5, log=T)
    l6 = dnbinom(bp2, size=size, mu=mu6, log=T)
    l7 = dnbinom(bh1, size=size, mu=mu7, log=T)
    l8 = dnbinom(bh2, size=size, mu=mu8, log=T)
    l = -sum(l1, l2, l3, l4, l5, l6, l7, l8)
    ifelse(is.nan(l), 100, l)
    #}}}
}
LL5 <- function(mu1,mu2,mu3,mu4,mu5,mu6,mu7,mu8) {
    #{{{
    l1 = dnbinom(p1, size=size, mu=mu1, log=T)
    l2 = dnbinom(p2, size=size, mu=mu2, log=T)
    l3 = dnbinom(h1, size=size, mu=mu3, log=T)
    l4 = dnbinom(h2, size=size, mu=mu4, log=T)
    l5 = dnbinom(bp1, size=size, mu=mu5, log=T)
    l6 = dnbinom(bp2, size=size, mu=mu6, log=T)
    l7 = dnbinom(bh1, size=size, mu=mu7, log=T)
    l8 = dnbinom(bh2, size=size, mu=mu8, log=T)
    l = -sum(l1, l2, l3, l4, l5, l6, l7, l8)
    ifelse(is.nan(l), 100, l)
    #}}}
}
    #}}}
    size = 1 / disp
    n_obs = length(p1)
    m1 = round(mean(p1)); m2 = round(mean(p2)); m3 = round(mean(h1)); m4 = round(mean(h2))
    m5 = round(mean(bp1)); m6 = round(mean(bp2)); m7 = round(mean(bh1)); m8 = round(mean(bh2))
    p1 = round(p1); p2 = round(p2); h1 = round(h1); h2 = round(h2)
    bp1 = round(bp1); bp2 = round(bp2); bh1 = round(bh1); bh2 = round(bh2)
    min_mu = 1e-2; max_mu = 1e8
    fit1 = mle2(LL1, start = list(mu1=m1,mu2=m2,mu3=m3,mu4=m4,mu5=m5,mu7=m7),
          lower = rep(min_mu, 6), upper=rep(max_mu, 6),
          method = method)#, nobs = n_obs)
    fit2 = mle2(LL2, start = list(mu1=m1,mu2=m2,mu3=m3,mu4=m4,mu5=m5,mu7=m7,mu8=m8),
          lower = rep(min_mu, 7), upper=rep(max_mu, 7),
          method = method)#, nobs = n_obs)
    fit3 = mle2(LL3, start = list(mu1=m1,mu2=m2,mu3=m3,mu4=m4,mu5=m5,mu6=m6,mu7=m7),
          lower = rep(min_mu, 7), upper=rep(max_mu, 7),
          method = method)#, nobs = n_obs)
    fit4 = mle2(LL4, start = list(mu1=m1,mu2=m2,mu3=m3,mu4=m4,mu5=m5,mu6=m6,mu7=m7),
          lower = rep(min_mu, 7), upper=rep(max_mu, 7),
          method = method)#, nobs = n_obs)
    fit5 = mle2(LL5, start = list(mu1=m1,mu2=m2,mu3=m3,mu4=m4,mu5=m5,mu6=m6,mu7=m7,mu8=m8),
          lower = rep(min_mu, 8), upper=rep(max_mu, 8),
          method = method)#, nobs = n_obs)
    #coef(fitc)
    bic = AICtab(fit1, fit2, fit3, fit4, fit5, k = log(n_obs), sort=F)
    tb = as_tibble(bic) %>%
        mutate(reg = c('conserved','unexpected','cis','trans','cis+trans')) %>%
        arrange(dAIC)
    tb$reg[[1]]
    #}}}
}
#}}}

i=98
#p1=ti$rc.p1[[i]]; p2=ti$rc.p2[[i]]; h1=ti$rc.h1[[i]]; h2=ti$rc.h2[[i]]; disp=ti$disp[i]

n_cpu = 8
#{{{ prepare for parallel computing
n_thread = n_cpu
cluster = new_cluster(n_thread)
cluster_library(cluster, "tidyverse")
cluster_library(cluster, "stats4")
cluster_library(cluster, "bbmle")
cluster_copy(cluster, 'calc_bic')
#}}}

#{{{ run cis/trans tests
fi = glue('{dirw}/01.raw.rds')
ra = readRDS(fi)
tmt = ra$tmt

min_rc = 10
ti = tmt %>%
    filter(trc.p1 + trc.p2 >= 2*min_rc, trc.h1+trc.h2 >= min_rc) %>%
    mutate(i= 1:n())

tw = ti %>%# dplyr::slice(1:10) %>%
    partition(cluster) %>%
    mutate(reg = pmap_chr(list(i,rc.p1,rc.p2,rc.h1,rc.h2,disp), calc_bic)) %>%
    collect()

to = tw %>% mutate(n.p1 = map_int(rc.p1, length), n.p2=map_int(rc.p2, length),
        n.h1 = map_int(rc.h1, length), n.h2=map_int(rc.h2, length)) %>%
    mutate(mrc.p1 = trc.p1/n.p1, mrc.p2 = trc.p2/n.p2,
           mrc.h1 = trc.h1/n.h1, mrc.h2 = trc.h2/n.h2) %>%
    mutate(prop.p = mrc.p1/(mrc.p1+mrc.p2), prop.h = mrc.h1/(mrc.h1+mrc.h2)) %>%
    select(cond,cross,gid,mrc.p1,mrc.p2,mrc.h1,mrc.h2,prop.p,prop.h,reg)

fo = file.path(dirw, '05.modes.rds')
saveRDS(to, fo)
#}}}



