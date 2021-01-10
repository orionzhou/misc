source('functions.R')
dirw = glue('{dird}/15_de')
drcs = c('up', 'down')
fmax=10

fi = glue("{dirw}/01.rds")
r = readRDS(fi)
th = r$th; tm = r$tm

#{{{ call DEGs
th1 = th %>% distinct(SampleID,tis,gt,cond)
cmp = th %>% distinct(tis) %>%
    mutate(cond1 = glue("{tis}.A632")) %>%
    mutate(cond2 = glue("{tis}.B73"))

x = call_deg(th1, tm, cmp)
fo = glue('{dirw}/01.de.rds')
saveRDS(x, fo)
#}}}



