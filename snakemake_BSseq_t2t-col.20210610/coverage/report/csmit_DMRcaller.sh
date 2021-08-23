#!/bin/bash

csmit -m 100G -c 47 "bash ./DMRcaller.sh 'WT_BSseq_Rep2_2013,WT_BSseq_Rep3_2013,WT_BSseq_Rep1_2014' 'cmt3_BSseq_Rep1' t2t-col.20210610 'Chr1,Chr2,Chr3,Chr4,Chr5' genomewide 6 CHG" & sleep 10;
csmit -m 100G -c 47 "bash ./DMRcaller.sh 'WT_BSseq_Rep2_2013,WT_BSseq_Rep3_2013,WT_BSseq_Rep1_2014' 'cmt3_BSseq_Rep1' t2t-col.20210610 'Chr1,Chr2,Chr3,Chr4,Chr5' genomewide 6 CHH" & sleep 10;

csmit -m 100G -c 25 "bash ./DMRcaller.sh 'WT_BSseq_Rep2_2013,WT_BSseq_Rep3_2013,WT_BSseq_Rep1_2014' 'kss_BSseq_Rep1' t2t-col.20210610 'Chr1,Chr2,Chr3,Chr4,Chr5' genomewide 6 CHG" & sleep 10;
csmit -m 100G -c 25 "bash ./DMRcaller.sh 'WT_BSseq_Rep2_2013,WT_BSseq_Rep3_2013,WT_BSseq_Rep1_2014' 'kss_BSseq_Rep1' t2t-col.20210610 'Chr1,Chr2,Chr3,Chr4,Chr5' genomewide 6 CHH" & sleep 10;

csmit -m 100G -c 25 "bash ./DMRcaller.sh 'WT_BSseq_Rep2_2013,WT_BSseq_Rep3_2013,WT_BSseq_Rep1_2014' 'cmt3_BSseq_Rep1' t2t-col.20210610 'Chr1,Chr2,Chr3,Chr4,Chr5' genomewide 6 CpG" & sleep 10;

csmit -m 100G -c 25 "bash ./DMRcaller.sh 'WT_BSseq_Rep2_2013,WT_BSseq_Rep3_2013,WT_BSseq_Rep1_2014' 'kss_BSseq_Rep1' t2t-col.20210610 'Chr1,Chr2,Chr3,Chr4,Chr5' genomewide 6 CpG" & sleep 10;
