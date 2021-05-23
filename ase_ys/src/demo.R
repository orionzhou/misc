require(tidyverse)
require(stringr)
require(glue)
dirp = "~/projects/misc/ase_ys"
dirw = glue('{dirp}/data')
setwd(dirw)
regs = c('conserved','unexpected','cis','trans','cis+trans')

tis = 'pericarp'; pre = 't1'
#{{{
fi = glue("{dirw}/{tis}.tsv")
t_rc = read_tsv(fi) %>%
    select(gid=GeneID, parA_rep1=B73.Pericarp.RawCounts.R1,
        parA_rep2=B73.Pericarp.RawCounts.R2,
        parA_rep3=B73.Pericarp.RawCounts.R3,
        parB_rep1=A632.Pericarp.RawCounts.R1,
        parB_rep2=A632.Pericarp.RawCounts.R2,
        parB_rep3=A632.Pericarp.RawCounts.R3,
        hybA_rep1=B73.allele.counts.R1,
        hybA_rep2=B73.allele.counts.R2,
        hybA_rep3=B73.allele.counts.R1_1,
        hybB_rep1=A632.allele.counts.R1,
        hybB_rep2=A632.allele.counts.R2,
        hybB_rep3=A632.allele.counts.R2_1)
#}}}

tis = 'silk'; pre = 't2'
#{{{
fi = glue("{dirw}/{tis}.tsv")
t_rc = read_tsv(fi) %>%
    select(gid=GeneID, parA_rep1=B73.silk.RawCounts.R1,
        parA_rep2=B73.silk.RawCounts.R2,
        parA_rep3=B73.silk.RawCounts.R3,
        parB_rep1=A632.silk.RawCounts.R1,
        parB_rep2=A632.silk.RawCounts.R2,
        parB_rep3=A632.silk.RawCounts.R3,
        hybA_rep1=B73.allele.counts.R1,
        hybA_rep2=B73.allele.counts.R2,
        hybA_rep3=B73.allele.counts.R1_1,
        hybB_rep1=A632.allele.counts.R1,
        hybB_rep2=A632.allele.counts.R2,
        hybB_rep3=A632.allele.counts.R2_1)
#}}}

t_sf = tibble(sid=colnames(t_rc)[-1]) %>%
    filter(!str_detect(sid, '^hybB')) %>%
    mutate(sid = str_replace(sid, 'hybA', 'hyb')) %>%
    mutate(sizeFactor=1)
t_dsp = t_rc %>% mutate(disp=0) %>% select(gid, disp)
#
fo = glue("{dirw}/{pre}_rc.tsv")
write_tsv(t_rc, fo)
fo = glue("{dirw}/{pre}_sf.tsv")
write_tsv(t_sf, fo)
fo = glue("{dirw}/{pre}_dsp.tsv")
write_tsv(t_dsp, fo)
system(glue("$git/demo/cis_trans/cis_trans.R --n_cpu 8 {pre}_rc.tsv {pre}_sf.tsv {pre}_dsp.tsv {pre}.tsv"))

