#!/bin/bash

# Usage on hydrogen node7:
# csmit -m 100G -c 47 "bash ./DMRcaller_hypo.sh 'WT_BSseq_Rep2_2013,WT_BSseq_Rep3_2013,WT_BSseq_Rep1_2014' 'cmt3_BSseq_Rep1' 't2t-col.20210610' 'Chr1,Chr2,Chr3,Chr4,Chr5' genomewide 6 CHG"

condition1=$1
condition2=$2
refbase=$3
chrName=$4
genomeRegion=$5
quantiles=$6
context=$7

source activate python

./DMRcaller_hypo.R --condition1 ${condition1} \
                   --condition2 ${condition2} \
                   --refbase ${refbase} \
                   --chrName ${chrName} \
                   --genomeRegion ${genomeRegion} \
                   --quantiles ${quantiles} \
                   --context ${context}

source deactivate
