#{{{ load & read
require(devtools)
require(GenomicFeatures)
load_all('~/git/rmaize')
require(progress)
require(ape)
require(ggtree)
require(ggforce)
require(Rtsne)
require(ggpubr)
require(lubridate)
options(dplyr.summarise.inform = F)
dirg = '~/data/genome'
dirp = "~/projects/misc/ase_ys"
dird = file.path(dirp, 'data')
dirf = file.path(dird, '95_figures', 'plots')
gcfg = read_genome_conf()
cols100 = colorRampPalette(rev(brewer.pal(n = 6, name = "RdYlBu")))(100)
cols100v = viridis_pal(direction=-1,option='magma')(100)
colbright <- function(col) {x = col2rgb(col); as.integer(.2126*x[1] + .7152*x[2] + .0722*x[3])}
cols36 = c(pal_ucscgb()(18)[8], pal_igv()(18), pal_ucscgb()(18)[c(1:7,9:18)])
brights36 = tibble(col=cols36) %>% mutate(bright=map_int(col,colbright)) %>%
    mutate(b = ifelse(bright<128, 'white','black')) %>% pull(b)
#}}}
gts3 = c("B73",'A632','B73xA632')
gt_map = list('B'='B73','A'='A632','BxA'='B73xA632')
tiss = c("glume",'pericarp','silk')

get_ds <- function(cond, condB, dds, gids) {
    #{{{
    res1 = results(dds, contrast = c("cond",cond,condB), pAdjustMethod='fdr')
    stopifnot(rownames(res1) == gids)
    tibble(gid = gids, padj = res1$padj, log2fc = res1$log2FoldChange) %>%
        replace_na(list(padj = 1))
    #}}}
}
call_deg_1 <- function(th, tm, base.cond = 'wt') {
    #{{{ th with 'cond1' and 'cond2' cols
    th1 = th %>% mutate(cond = str_c(cond1, cond2, sep='.'))
    tm1 = tm %>% filter(SampleID %in% th1$SampleID)
    ct = th1 %>% distinct(cond, cond1, cond2) %>% filter(cond2 != base.cond) %>%
        mutate(condB = str_c(cond1, base.cond, sep="."))
    #{{{ prepare data
    vh = th1 %>% mutate(cond2 = factor(cond2)) %>% arrange(SampleID)
    vh.d = column_to_rownames(as.data.frame(vh), var = 'SampleID')
    gids = tm1 %>% group_by(gid) %>% summarise(n.sam = sum(ReadCount >= 10)) %>%
        filter(n.sam > .2 * nrow(vh)) %>% pull(gid)
    vm = tm1 %>% filter(gid %in% gids) %>%
        select(SampleID, gid, ReadCount)
    x = readcount_norm(vm)
    mean.lib.size = mean(x$tl$libSize)
    vm = x$tm
    vm.w = vm %>% select(SampleID, gid, ReadCount) %>% spread(SampleID, ReadCount)
    vm.d = column_to_rownames(as.data.frame(vm.w), var = 'gid')
    stopifnot(identical(rownames(vh.d), colnames(vm.d)))
    #}}}
    # DESeq2
    dds = DESeqDataSetFromMatrix(countData=vm.d, colData=vh.d, design=~cond)
    dds = estimateSizeFactors(dds)
    dds = estimateDispersions(dds, fitType = 'parametric')
    disp = dispersions(dds)
    #dds = nbinomLRT(dds, reduced = ~ 1)
    dds = nbinomWaldTest(dds)
    resultsNames(dds)
    res = ct %>% mutate(ds = map2(cond, condB, get_ds, dds = dds, gids = gids))
    res %>% select(cond1, cond, condB, ds)
    #}}}
}
call_deg <- function(th, tm, comps) {
    #{{{ th with 'cond' column and comps with 'cond1' and 'cond2' columns
    th1 = th
    tm1 = tm %>% filter(SampleID %in% th1$SampleID)
    #{{{ prepare data
    vh = th1 %>% arrange(SampleID)
    vh.d = column_to_rownames(as.data.frame(vh), var = 'SampleID')
    gids = tm1 %>% group_by(gid) %>% summarise(n.sam = sum(ReadCount >= 10)) %>%
        filter(n.sam > .2 * nrow(vh)) %>% pull(gid)
    vm = tm1 %>% filter(gid %in% gids) %>%
        select(SampleID, gid, ReadCount)
    x = readcount_norm(vm)
    mean.lib.size = mean(x$tl$libSize)
    vm = x$tm
    vm.w = vm %>% select(SampleID, gid, ReadCount) %>% spread(SampleID, ReadCount)
    vm.d = column_to_rownames(as.data.frame(vm.w), var = 'gid')
    stopifnot(identical(rownames(vh.d), colnames(vm.d)))
    #}}}
    # DESeq2
    dds = DESeqDataSetFromMatrix(countData=vm.d, colData=vh.d, design=~cond)
    dds = estimateSizeFactors(dds)
    dds = estimateDispersions(dds, fitType = 'parametric')
    disp = dispersions(dds)
    #dds = nbinomLRT(dds, reduced = ~ 1)
    dds = nbinomWaldTest(dds)
    resultsNames(dds)
    res = comps %>% mutate(ds = map2(cond1, cond2, get_ds, dds=dds, gids=gids))
    res
    #}}}
}



