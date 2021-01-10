source('functions.R')
dirw = glue('{dird}/10_qc')

#{{{ read in & format
yid = 'rn20j'
res = rnaseq_cpm(yid)
#
th = res$th %>% mutate(gt=str_replace(Genotype,'A632d','A632')) %>%
    mutate(tis = factor(Tissue, levels=tiss)) %>%
    mutate(gt = factor(gt, levels=gts3)) %>%
    arrange(tis, gt, SampleID) %>%
    mutate(cond=glue("{tis}.{gt}")) %>%
    mutate(cond = as_factor(cond)) %>%
    mutate(SampleID = as_factor(SampleID)) %>%
    select(SampleID,tis,gt,cond)
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(SampleID = factor(SampleID, levels=levels(th$SampleID))) %>%
    mutate(value=asinh(CPM))

fo = glue("{dirw}/01.rds")
r = list(th = th, tm = tm)
saveRDS(r, fo)
#}}}

