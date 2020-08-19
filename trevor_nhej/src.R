require(tidyverse)
dirw = '/home/springer/zhoux379/projects/misc/trevor_nhej'
#{{{
get_longest_del <- function(s)
    str_locate_all(s, "[ATCGN]-+[ATCGN]")[[1]] %>% as_tibble() %>%
        mutate(size=end-start-1) %>% arrange(-size) %>% slice(1)
max_match_nt <- function(s1, s2, direction=1) {
    #{{{
    len1 = nchar(s1); len2 = nchar(s2)
    len = min(len1, len2)
    cnt = 0
    for (i in 1:len) {
        j1 = i; j2 = i
        if (direction == -1) {
            j1 = len1 - i + 1;
            j2 = len2 - i + 1
        }
        if (str_sub(s1, j1, j1) != str_sub(s2, j2, j2)) { break }
        cnt = cnt + 1
    }
    cnt
    #}}}
}
run_mmej_pipe <- function(pre, diri, diro) {
    #{{{
    if(!dir.exists(diro)) system(sprintf("mkdir -p %s", diro))
    fi = sprintf("%s/%s.txt", diri, pre)
    ti = read_tsv(fi) %>%
        rename(aSeq=1,rSeq=2,ref=3,status=4,n_del=5,n_ins=6,n_mut=7,num_reads=8,pct_reads=9)
    ti2 = ti %>%
        mutate(dpos = map(aSeq, get_longest_del))
    #
    ti3 = ti2 %>% unnest(dpos) %>%
        mutate(dseq = str_sub(rSeq, start+1, end-1)) %>%
        mutate(lstart=start-19, lend=start, rstart=end, rend=end+19) %>%
        mutate(lseq = str_sub(aSeq, lstart, lend)) %>%
        mutate(rseq = str_sub(aSeq, rstart, rend)) %>%
        mutate(lnt = map2_dbl(dseq, lseq, max_match_nt, direction=-1)) %>%
        mutate(rnt = map2_dbl(dseq, rseq, max_match_nt, direction=1))
    #
    to = ti3 %>% mutate(l_bp_hom = lnt, r_bp_hom = rnt, bp_hom = pmax(lnt, rnt)) %>%
        mutate(type=ifelse(bp_hom<=1, 'NHEJ', 'MMEJ')) %>%
        select(-start,-end,-lstart,-lend,-rstart,-rend,-lseq,-rseq,-lnt,-rnt)
    #
    fo = sprintf("%s/%s.txt", diro, pre)
    write_tsv(to, fo)
    #}}}
}
#}}}

#{{{ set 1
diri = file.path(dirw, 'MMEJ')
diro = file.path(dirw, 'MMEJ_out')
pres = c(
    "pMG198_ms26_gRNA1",
    "pMG198_ms45_gRNA1",
    "pMG199_ms26_gRNA2",
    "pMG199_ms45_gRNA2",
    "pMG201_ms26_gRNA1",
    "pMG202_ms26_gRNA2",
    "pMG202_ms45_gRNA2"
)
#}}}

#{{{ set 2
diri = file.path(dirw, 'mmej_020720')
diro = file.path(dirw, 'mmej_020720_out')
pres = c(
    "Ms26_gRNA1_cas9",
    "Ms26_gRNA1_trex2_cas9",
    "Ms26_gRNA2_cas9_copy",
    "Ms26_gRNA2_cas9",
    "Ms45_gRNA1_cas9",
    "Ms45_gRNA1_trex2_cas9",
    "Ms45_gRNA2_trex2_cas9"
)
#}}}

#{{{ set3
diri = file.path(dirw, 'mmej_032520')
diro = file.path(dirw, 'mmej_032520_out')
pres = sprintf("5907-%d", c(1:5,7:9))
#}}}

tibble(pre = pres) %>% mutate(x = map(pre, run_mmej_pipe, diri=diri, diro=diro))
