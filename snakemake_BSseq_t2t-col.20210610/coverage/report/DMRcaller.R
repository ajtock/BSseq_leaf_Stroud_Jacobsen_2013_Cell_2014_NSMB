#!/usr/bin/env Rscript

# Use DMRcaller v1.24.0 to identify DMRs between two conditions
# (e.g. mutant and wild type)


# Usage:
# conda activate python
# ./DMRcaller.R --condition1 'WT_BSseq_Rep1_2014,WT_BSseq_Rep2_2013,WT_BSseq_Rep3_2013' \
#               --condition2 'cmt3_BSseq_Rep1' \
#               --refbase 't2t-col.20210610' \
#               --chrName 'Chr1,Chr2,Chr3,Chr4,Chr5' \
#               --context 'CHG'
# conda deactivate

library(argparse)
library(DMRcaller)
library(betareg) # required by DMRcaller::computeDMRsReplicates 

# Create parser object
parser <- ArgumentParser()

# Specify arguments
# ArgumentParser will add a help option by default
parser$add_argument("--condition1", type = "character", default = "WT_BSseq_Rep1_2014,WT_BSseq_Rep2_2013,WT_BSseq_Rep3_2013",
                    help="Sample replicate prefixes for first condition for inclusion in DMRcaller contrast. Default: %(default)s")
parser$add_argument("--condition2", type = "character",
                    help="Sample replicate prefixes for second condition for inclusion in DMRcaller contrast.")
parser$add_argument("--refbase", type = "character", default = "t2t-col.20210610",
                    help="Reference genome base name. Default: %(default)s")
parser$add_argument("--chrName", type = "character", default = "Chr1,Chr2,Chr3,Chr4,Chr5",
                    help="Reference genome chromosome names for inclusion in DMRcaller contrast. Default: %(default)s")
parser$add_argument("--context", type = "character",
                    help="cytosine methylation context.")

# parse arguments
args <- parser$parse_args()
args_file <- "tempArgsObjectFile.rds"
#saveRDS(args, args_file); print(args)

#system("./DMRcaller.R --condition1 'WT_BSseq_Rep1_2014,WT_BSseq_Rep2_2013,WT_BSseq_Rep3_2013' --condition2 'cmt3_BSseq_Rep1' --refbase t2t-col.20210610 --chrName 'Chr1,Chr2,Chr3,Chr4,Chr5' --context CHG")

args <- readRDS(args_file)
args$condition1 <- unlist(strsplit(args$condition1, split = ","))
args$condition2 <- unlist(strsplit(args$condition2, split = ","))
args$chrName <- unlist(strsplit(args$chrName, split = ","))

condition1_Reps <- mclapply(seq_along(args$condition1), function(x) {
  readBismark(paste0(args$condition1[x], "_MappedOn_", args$refbase, "_dedup_", args$context, ".CX_report.txt.gz"))
}, mc.cores = length(args$condition1))
condition2_Reps <- mclapply(seq_along(args$condition2), function(x) {
  readBismark(paste0(args$condition2[x], "_MappedOn_", args$refbase, "_dedup_", args$context, ".CX_report.txt.gz"))
}, mc.cores = length(args$condition2))

## For pooled-replicate analysis
#condition1_Reps_pooled <- poolMethylationDatasets(GRangesList(condition1_Reps))
#condition2_Reps_pooled <- poolMethylationDatasets(GRangesList(condition2_Reps))

# For biological replicate analysis
# (However, computeDMRsReplicates fails with error message:
# "Error in apply(readsM2, 1, sum) : dim(X) must have a positive length")
conditions_Reps_list <- c(condition1_Reps, condition2_Reps)

# Combine replicates into a single GRanges object containing
# data for each condition and replicate
joined_Reps <- joinReplicates(methylationData1 = conditions_Reps_list[[1]],
                              methylationData2 = conditions_Reps_list[[2]],
                              usecomplete = FALSE)
for(x in 3:length(conditions_Reps_list)) {
  joined_Reps <- joinReplicates(methylationData1 = joined_Reps,
                                methylationData2 = conditions_Reps_list[[x]]) 
}

# Get ranges corresponding to the given context
joined_Reps <- joined_Reps[joined_Reps$context == sub("p", "", args$context)]

# Get ranges corresponding to those in chrName
#joined_Reps <- joined_Reps[seqnames(joined_Reps) %in% args$chrName]
joined_Reps <- keepSeqlevels(joined_Reps, args$chrName, pruning.mode="coarse")

# Sort by seqnames, strand, start, and end
joined_Reps <- sortSeqlevels(joined_Reps)
joined_Reps <- sort(joined_Reps, ignore.strand = TRUE)

pdf(paste0("plotMethylationDataCoverage_",
           args$condition1[1], "_", args$condition1[2],
           args$condition1[3], "_", args$condition2[1],
           ".pdf"))
plotMethylationDataCoverage(methylationData1 = conditions_Reps_list[[1]],
                            methylationData2 = conditions_Reps_list[[2]],
                            breaks = c(1, 3, 6, 9, 12, 15),
                            regions = NULL,
                            conditionsNames = c(args$condition1[1], args$condition1[2]),
                            context = sub("p", "", args$context),
                            proportion = TRUE,
                            labels=LETTERS,
                            contextPerRow = FALSE)
plotMethylationDataCoverage(methylationData1 = conditions_Reps_list[[3]],
                            methylationData2 = conditions_Reps_list[[4]],
                            breaks = c(1, 3, 6, 9, 12, 15),
                            regions = NULL,
                            conditionsNames = c(args$condition1[3], args$condition2[1]),
                            context = sub("p", "", args$context),
                            proportion = TRUE,
                            labels=LETTERS,
                            contextPerRow = FALSE)
dev.off()


# Create condition vector
joined_Reps_conditions <- gsub("_.+", "", c(args$condition1, args$condition2))

print("joined_Reps conditions:")
print(joined_Reps_conditions)

if(args$context == "CpG") {
  minProportionDifference_context <- 0.4
} else if(args$context == "CHG") {
  minProportionDifference_context <- 0.2
} else if(args$context == "CHH") {
  minProportionDifference_context <- 0.1
}
print(paste0(args$context, " minProportionDifference = ", minProportionDifference_context))

# Compute DMRs using "bins" method
DMRs_bins <- computeDMRs(methylationData1 = ,
                         methylationData2 = ,
                         regions = NULL,
                         context = sub("p", "", args$context),
                         method = "bins",
                         binSize = 100,
                         test = "fisher",
                         pValueThreshold = 0.01,
                         minCytosinesCount = 4,
                         minProportionDifference = minProportionDifference_context,
                         minGap = 200,
                         minSize = 50,
                         minReadsPerCytosine = 4,
                         cores = 1)

# Compute DMRs using "bins" method
DMRsReplicates_bins <- computeDMRsReplicates(methylationData = joined_Reps,
                                             condition = joined_Reps_conditions,
                                             regions = NULL,
                                             context = sub("p", "", args$context),
                                             method = "bins",
                                             binSize = 100,
                                             test = "betareg",
                                             pseudocountM = 1,
                                             pseudocountN = 2,
                                             pValueThreshold = 0.01,
                                             minCytosinesCount = 4,
                                             minProportionDifference = minProportionDifference_context,
                                             minGap = 200,
                                             minSize = 50,
                                             minReadsPerCytosine = 4,
                                             cores = 1)

DMRsReplicates_bins <- computeDMRsReplicates(methylationData = methylationData,
                                             condition = c("WT", "WT", "mut", "mut"),
                                             regions = NULL,
                                             context = "CHH",
                                             method = "bins",
                                             binSize = 100,
                                             test = "betareg",
                                             pseudocountM = 1,
                                             pseudocountN = 2,
                                             pValueThreshold = 0.01,
                                             minCytosinesCount = 4,
                                             minProportionDifference = 0.1,
                                             minGap = 200,
                                             minSize = 50,
                                             minReadsPerCytosine = 4,
                                             cores = 1)

## Compute DMRs using "neighbourhood" method
#DMRsReplicates_neighbourhood <- computeDMRsReplicates(methylationData = joined_Reps,
#                                                      condition = joined_Reps_conditions,
#                                                      regions = NULL,
#                                                      context = sub("p", "", args$context),
#                                                      method = "neighbourhood",
#                                                      test = "betareg",
#                                                      pseudocountM = 1,
#                                                      pseudocountN = 2,
#                                                      pValueThreshold = 0.01,
#                                                      minCytosinesCount = 4,
#                                                      minProportionDifference = minProportionDifference_context,
#                                                      minGap = 200,
#                                                      minSize = 50,
#                                                      minReadsPerCytosine = 4,
#                                                      cores = detectCores())

