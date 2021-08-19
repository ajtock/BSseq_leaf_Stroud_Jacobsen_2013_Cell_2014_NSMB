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
library(rtracklayer)

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
parser$add_argument("--chrName", type = "character",
                    help="Reference genome chromosome names for inclusion in DMRcaller contrast.")
parser$add_argument("--context", type = "character",
                    help="cytosine methylation context.")

# parse arguments
args <- parser$parse_args()
args_file <- "tempArgsObjectFile.rds"
#saveRDS(args, args_file); print(args)

#system("./DMRcaller.R --condition1 'WT_BSseq_Rep1_2014,WT_BSseq_Rep2_2013,WT_BSseq_Rep3_2013' --condition2 'cmt3_BSseq_Rep1' --refbase t2t-col.20210610 --chrName 'Chr4' --context CHG")
#system("./DMRcaller.R --condition1 'WT_BSseq_Rep1_2014,WT_BSseq_Rep2_2013,WT_BSseq_Rep3_2013' --condition2 'cmt3_BSseq_Rep1' --refbase t2t-col.20210610 --chrName 'Chr1,Chr2,Chr3,Chr4,Chr5' --context CHG")

args <- readRDS(args_file)
args$condition1 <- unlist(strsplit(args$condition1, split = ","))
args$condition2 <- unlist(strsplit(args$condition2, split = ","))
args$chrName <- unlist(strsplit(args$chrName, split = ","))

# Set context-specific minProportionDifference to pass to computeDMRs or computeDMRsReplicates
if(args$context == "CpG") {
  minProportionDifference_context <- 0.4
} else if(args$context == "CHG") {
  minProportionDifference_context <- 0.2
} else if(args$context == "CHH") {
  minProportionDifference_context <- 0.1
}
print(paste0(args$context, " minProportionDifference = ", minProportionDifference_context))

# Load methylation data
condition1_Reps <- mclapply(seq_along(args$condition1), function(x) {
  readBismark(paste0(args$condition1[x], "_MappedOn_", args$refbase, "_dedup_", args$context, ".CX_report.txt.gz"))
}, mc.cores = length(args$condition1))
condition2_Reps <- mclapply(seq_along(args$condition2), function(x) {
  readBismark(paste0(args$condition2[x], "_MappedOn_", args$refbase, "_dedup_", args$context, ".CX_report.txt.gz"))
}, mc.cores = length(args$condition2))

# Calculate conversion rate and adjust coverage values accordingly
# See https://doi.org/10.1007/978-1-0716-1134-0_21
# (full link: https://link.springer.com/protocol/10.1007%2F978-1-0716-1134-0_21 )
for(x in 1:length(condition1_Reps)) {
  print(paste0(x, ">>"))
  # Get ranges corresponding to the given context
  condition1_Reps[[x]] <- condition1_Reps[[x]][condition1_Reps[[x]]$context == sub("p", "", args$context)]

  # Extract chloroplast methylation data
  condition1_Repx_ChrC <- condition1_Reps[[x]][condition1_Reps[[x]]$context == sub("p", "", args$context) &
                                               seqnames(condition1_Reps[[x]]) == "ChrC"]

  # Calculate conversion
  condition1_Repx_conv <- 1 - ( sum(mcols(condition1_Repx_ChrC)$readsM) / sum(mcols(condition1_Repx_ChrC)$readsN) )

  # Adjust methylated and total read counts (readsM and readsN)
  condition1_Reps[[x]]$readsM <- round( condition1_Reps[[x]]$readsM - condition1_Reps[[x]]$readsN * (1 - condition1_Repx_conv) )
  condition1_Reps[[x]]$readsM[condition1_Reps[[x]]$readsM < 0] <- 0 
  condition1_Reps[[x]]$readsN <- round( condition1_Reps[[x]]$readsN * condition1_Repx_conv )

  # Remove superfluous ranges
  # Get ranges corresponding to those in chrName
  condition1_Reps[[x]] <- keepSeqlevels(condition1_Reps[[x]], args$chrName, pruning.mode = "coarse")

  # Sort by seqnames, start and end
  condition1_Reps[[x]] <- sortSeqlevels(condition1_Reps[[x]])
  condition1_Reps[[x]] <- sort(condition1_Reps[[x]], ignore.strand = TRUE)

  print(condition1_Reps[[x]])
  print(paste0("<<", x))
}

for(x in 1:length(condition2_Reps)) {
  print(paste0(x, ">>"))
  # Get ranges corresponding to the given context
  condition2_Reps[[x]] <- condition2_Reps[[x]][condition2_Reps[[x]]$context == sub("p", "", args$context)]

  # Extract chloroplast methylation data
  condition2_Repx_ChrC <- condition2_Reps[[x]][condition2_Reps[[x]]$context == sub("p", "", args$context) &
                                               seqnames(condition2_Reps[[x]]) == "ChrC"]

  # Calculate conversion
  condition2_Repx_conv <- 1 - ( sum(mcols(condition2_Repx_ChrC)$readsM) / sum(mcols(condition2_Repx_ChrC)$readsN) )

  # Adjust methylated and total read counts (readsM and readsN)
  condition2_Reps[[x]]$readsM <- round( condition2_Reps[[x]]$readsM - condition2_Reps[[x]]$readsN * (1 - condition2_Repx_conv) )
  condition2_Reps[[x]]$readsM[condition2_Reps[[x]]$readsM < 0] <- 0 
  condition2_Reps[[x]]$readsN <- round( condition2_Reps[[x]]$readsN * condition2_Repx_conv )

  # Remove superfluous ranges
  # Get ranges corresponding to those in chrName
  condition2_Reps[[x]] <- keepSeqlevels(condition2_Reps[[x]], args$chrName, pruning.mode = "coarse")

  # Sort by seqnames, start and end
  condition2_Reps[[x]] <- sortSeqlevels(condition2_Reps[[x]])
  condition2_Reps[[x]] <- sort(condition2_Reps[[x]], ignore.strand = TRUE)

  print(condition2_Reps[[x]])
  print(paste0("<<", x))
}


if(length(condition2_Reps) == 1) {

  # Plot methylation coverage for each replicate and condition
  pdf(paste0("plotMethylationDataCoverage_",
             paste0(args$condition1, collapse = "_"), "_",
             paste0(args$condition2, collapse = "_"), "_",
             paste0(args$chrName, collapse = "_"), "_",
             paste0(args$context, collapse = "_"),
             ".pdf"))
  plotMethylationDataCoverage(methylationData1 = condition1_Reps[[1]],
                              methylationData2 = condition1_Reps[[2]],
                              breaks = c(1, 3, 6, 9, 12, 15),
                              regions = NULL,
                              conditionsNames = c(args$condition1[1], args$condition1[2]),
                              context = sub("p", "", args$context),
                              proportion = TRUE,
                              labels=LETTERS,
                              contextPerRow = FALSE)
  plotMethylationDataCoverage(methylationData1 = condition1_Reps[[3]],
                              methylationData2 = condition2_Reps[[1]],
                              breaks = c(1, 3, 6, 9, 12, 15),
                              regions = NULL,
                              conditionsNames = c(args$condition1[3], args$condition2[1]),
                              context = sub("p", "", args$context),
                              proportion = TRUE,
                              labels=LETTERS,
                              contextPerRow = FALSE)
  dev.off()

  # Compute DMRs using "bins" method
  DMRs_per_Rep_list_bins <- lapply(seq_along(condition1_Reps), function(x) {
    computeDMRs(methylationData1 = condition1_Reps[[x]],
                methylationData2 = condition2_Reps[[1]] ,
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
                cores = 48) 
  })

  DMRs_Rep1_bins_DMRs_Repx_bins_overlap <- lapply(seq_along(length(DMRs_per_Rep_list_bins), function(x) {
    findOverlaps(query = DMRs_per_Rep_list_bins[[1]],
                 subject = DMRs_per_Rep_list_bins[[x]],
                 type = "any", select = "all",
                 ignore.strand = FALSE)
  })

stopifnot( length( findOverlaps(query = ranLocAcc1GR,
                                subject = genomeMaskGR,
                                type = "any", select = "all",
                                ignore.strand = TRUE) ) == 0 )

    

} else {
  # For biological replicate analysis using computeDMRsReplicates:
  # (Note that this requires equal numbers of replicates for each condition,
  # or > 1 replicate for each condition. If this is not the case,
  # computeDMRsReplicates fails with an error message:
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
  
  ## Get ranges corresponding to the given context
  #joined_Reps <- joined_Reps[joined_Reps$context == sub("p", "", args$context)]
  #
  ## Get ranges corresponding to those in chrName
  ##joined_Reps <- joined_Reps[seqnames(joined_Reps) %in% args$chrName]
  #joined_Reps <- keepSeqlevels(joined_Reps, args$chrName, pruning.mode="coarse")
  #
  ## Sort by seqnames, start and end
  #joined_Reps <- sortSeqlevels(joined_Reps)
  #joined_Reps <- sort(joined_Reps, ignore.strand = TRUE)
  #
  #print(joined_Reps)
  
  
  # Create condition vector
  joined_Reps_conditions <- gsub("_.+", "", c(args$condition1, args$condition2))
  
  print("joined_Reps conditions:")
  print(joined_Reps_conditions)
  
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
  
}


## For pooled-replicate analysis
#condition1_Reps_pooled <- poolMethylationDatasets(GRangesList(condition1_Reps))
#condition2_Reps_pooled <- poolMethylationDatasets(GRangesList(condition2_Reps))
