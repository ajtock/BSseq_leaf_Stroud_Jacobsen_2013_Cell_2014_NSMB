#!/applications/R/R-4.0.0/bin/Rscript

# Use DMRcaller v1.20.0 to identify DMRs between two conditions
# (e.g. mutant and wild type)


# Usage:
# ./DMRcaller.R --condition1 'WT_BSseq_Rep1_2014,WT_BSseq_Rep2_2013,WT_BSseq_Rep3_2013' \
#               --condition2 'met1_BSseq_Rep1' \
#               --refbase 't2t-col.20210610' \
#               --context 'CHG'

#library(DMRcaller)
library(argparse)

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
parser$add_argument("--context", type = "character",
                    help="cytosine methylation context.")

# parse arguments
args <- parser$parse_args()
args_file <- "tempArgsObjectFile.rds"
saveRDS(args, args_file); print(args)

#system("./DMRcaller.R --condition1 'WT_BSseq_Rep1_2014,WT_BSseq_Rep2_2013,WT_BSseq_Rep3_2013' --condition2 'met1_BSseq_Rep1' --refbase t2t-col.20210610 --context CHG")

args <- readRDS(args_file)
args$condition1 <- unlist(strsplit(args$condition1, split = ","))
args$condition2 <- unlist(strsplit(args$condition2, split = ","))


#mCon1_Rep1 <- readBismark("WT_BSseq_Rep1_2014_MappedOn_t2t-col.20210610_dedup_CHH.CX_report.txt.gz")
#mCon1_Rep2 <- readBismark("WT_BSseq_Rep2_2013_MappedOn_t2t-col.20210610_dedup_CHH.CX_report.txt.gz")
#mCon1_Rep3 <- readBismark("WT_BSseq_Rep3_2013_MappedOn_t2t-col.20210610_dedup_CHH.CX_report.txt.gz")



condition1_Reps <- mclapply(seq_along(args$condition1), function(x) {
  readBismark(paste0(args$condition1[x], "_MappedOn_", args$refbase, "_dedup_", args$context, ".CX_report.txt.gz"))
}, mc.cores = length(args$condition1))
condition2_Reps <- mclapply(seq_along(args$condition2), function(x) {
  readBismark(paste0(args$condition2[x], "_MappedOn_", args$refbase, "_dedup_", args$context, ".CX_report.txt.gz"))
}, mc.cores = length(args$condition2))

conditions_Reps_list <- c(condition1_Reps, condition2_Reps)

# Combine replicates into a single GRanges object containing
# data for each condition and replicate
repsJoined <- joinReplicates(methylationData1 = conditions_Reps_list[[1]],
                             methylationData2 = conditions_Reps_list[[2]],
                             usecomplete = FALSE)
for(x in 3:length(conditions_Reps_list)) {
  repsJoined <- joinReplicates(methylationData1 = repsJoined,
                               methylationData2 = conditions_Reps_list[[x]]) 
}

# Get ranges corresponding to the given context
repsJoined <- repsJoined[repsJoined$context == sub("p", "", args$context)]

# Sort by seqnames, strand, start, and end
repsJoined <- sortSeqlevels(repsJoined)
repsJoined <- sort(repsJoined)




