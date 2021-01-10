source('functions.R')
dirw = glue('{dird}/16_ase')
regs = c('conserved','unexpected','cis','trans','cis+trans')
crosses = c("B73xA632")

#{{{ organize data for tests
fi = glue("{dirw}/01.rds")
r = readRDS(fi)
th = r$th %>% mutate(cond=tis)
tm = r$tm

#{{{ get sizeFactor and dispersions using DESeq2
require(DESeq2)
tm = r$tm %>% select(gid, SampleID, ReadCount)
tw = tm %>%
    spread(SampleID, ReadCount) %>%
    replace(., is.na(.), 0)
gids = tw$gid
twd = data.frame(tw[,-1])
rownames(twd) = tw$gid
th2 = th %>% mutate(sid = SampleID, gt=factor(gt), cond=factor(cond))
thd = column_to_rownames(as.data.frame(th2), var = 'sid')
stopifnot(identical(rownames(thd), colnames(twd)))
dds = DESeqDataSetFromMatrix(countData=twd, colData=thd, design = ~ cond)
dds = estimateSizeFactors(dds)
sf = sizeFactors(dds)
t_sf = tibble(SampleID = names(sf), sizeFactor = as.numeric(sf))
dds = estimateDispersions(dds, fitType = 'parametric')
stopifnot(identical(gids, names(dds)))
t_disp = tibble(gid=gids, disp = dispersions(dds))
tl = th %>% inner_join(t_sf, by = 'SampleID')
#}}}

th1 = th %>%
    filter(str_detect(gt, 'x')) %>% distinct(cond, gt) %>%
    mutate(gt = as.character(gt)) %>% mutate(cross = gt) %>%
    dplyr::rename(h=gt) %>% separate(h, c("p1",'p2'), sep='x', remove=F) %>%
    gather(gen, gt, -cond, -cross) %>%
    inner_join(th, by=c("cond",'gt')) %>% arrange(cond,cross,gen,gt,SampleID)

tm_p = th1 %>% inner_join(t_sf, by='SampleID') %>%
    inner_join(tm, by=c('SampleID')) %>%
    mutate(nRC = ReadCount / sizeFactor) %>%
    group_by(cond,cross,gen,gid) %>%
    summarise(rc=list(nRC), trc=sum(nRC)) %>% ungroup() %>%
    pivot_wider(names_from=gen, values_from=c(rc,trc), names_sep='.')
#
tm_ase = th1 %>% filter(gen=='h') %>%
    inner_join(t_sf, by='SampleID') %>%
    inner_join(res$ase_gene, by=c("SampleID")) %>%
    mutate(nRC1 = allele1 / sizeFactor, nRC2 = allele2 / sizeFactor) %>%
    group_by(cond,cross,gid) %>%
    summarise(rc.h1=list(nRC1), rc.h2=list(nRC2),
              trc.h1=sum(nRC1), trc.h2=sum(nRC2)) %>% ungroup()
#
tmt = tm_p %>% inner_join(t_disp, by='gid') %>%
    inner_join(tm_ase, by=c("cond",'cross','gid')) %>%
    select(cond,cross,gid,disp,trc.p1,trc.p2,trc.h,trc.h1,trc.h2,
           rc.p1,rc.p2,rc.h,rc.h1,rc.h2)

ra = list(tl=tl, disp=t_disp, th=th1, th.cmp = th2, tmt = tmt)
fo = file.path(dirw, '01.raw.rds')
saveRDS(ra, fo)
#}}}

# run st.16.ase.1.R

#{{{ simple cis/trans
fi = file.path(dirw, '05.modes.rds')
ti = readRDS(fi)
ti %>% count(cond,cross,reg) %>% spread(reg,n) %>% print(n=40)
ct_basic = ti

nconds = crossing(x=tiss, y=crosses) %>% mutate(cond=str_c(x,y,sep='_')) %>%
    mutate(x=factor(x,levels=tiss), y=factor(y,levels=crosses)) %>%
    arrange(x,y) %>% pull(cond)
tpx = tibble(reg=regs) %>% mutate(x=1:n())
tp = ti %>% mutate(cond = str_c(cond, cross, sep='_')) %>%
    inner_join(tpx, by='reg') %>%
    mutate(cond = factor(cond, levels=nconds)) %>%
    mutate(reg=factor(reg, levels=regs))
linecol = 'azure3'; lty = 'solid'
cols9 = pal_npg()(10)
#{{{ scatter plot
pa = ggplot(tp, aes(x=prop.p, y=prop.h)) +
    #geom_vline(xintercept=0, linetype=lty, color=linecol) +
    #geom_hline(yintercept=0, linetype=lty, color=linecol) +
    geom_point(aes(color=reg, shape=reg), size=1) +
    geom_abline(intercept=0, slope=1, linetype=lty, color=linecol) +
    scale_x_continuous(name='Allele 1 proportion in Parent', limits=c(0,1),expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(name='Allele 1 proportion in F1', limits=c(0,1),expand=expansion(mult=c(.05,.05))) +
    scale_color_manual(name='Mode:', values=cols9) +
    scale_shape_manual(name='Mode:', values=c(1:5)) +
    facet_wrap(~cond, ncol=1) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
           legend.vjust = -.4, legend.box='h',
           xtick=T, xtext=T, xtitle=T, ytitle=T, ytick=T, ytext=T) +
    guides(color=guide_legend(nrow=2))
#}}}
fp = sprintf('%s/10.modes.pdf', dirw)
ggsave(pa, file = fp, width = 4, height = 8)

#{{{ bar plot showing counts
tp1 = tp %>% count(cond, x, reg)
tp1s = tp1 %>% group_by(cond) %>% summarise(n=sum(n)) %>% ungroup() %>%
    mutate(lab=sprintf("N=%d", n))
ymax = max(tp1$n) * 1.05
pb = ggplot(tp1) +
    geom_bar(aes(x=x, y=n, fill=reg), width=.7, stat='identity') +
    geom_text(data=tp1s,aes(x=5.6,y=ymax,label=lab), size=2.5, vjust=.5,hjust=1) +
    geom_text(aes(x=x, y=n+100, label=n), size=2.5, vjust=0) +
    scale_x_continuous(breaks=tpx$x, labels=tpx$reg, expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(name='Number Genes', expand=expansion(mult=c(.05,.05))) +
    scale_fill_manual(name='Modes:', values=cols9) +
    facet_wrap(~cond, ncol=1) +
    otheme(legend.pos='none', legend.dir='v', legend.title=T,
           xtick=T, xtext=T, xtitle=F, ytitle=T, ytick=T, ytext=T) +
    theme(axis.text.x = element_text(angle=15, hjust=1, vjust=1))
#}}}
fp = sprintf('%s/10.modes.cnt.pdf', dirw)
ggsave(pb, file = fp, width = 3, height = 8)

fo = glue("{dirw}/11.modes.pdf")
ggarrange(pa, pb, nrow=1, ncol=2, labels=LETTERS[1:3], widths=c(1,1)) %>%
    ggexport(filename=fo, width=6, height=8)
#}}}

fo = glue('{dirw}/20.rds')
saveRDS(ct_basic, fo)
fo = glue('{dirw}/20.tsv.gz')
write_tsv(ct_basic, fo)


