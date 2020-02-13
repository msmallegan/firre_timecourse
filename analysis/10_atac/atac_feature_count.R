#!/usr/bin/env Rscript

knitr::opts_chunk$set(echo = TRUE)
options(stringsAsFactors = FALSE)
library(tidyverse)
library(Rsubread)

# let's grab all the bam files. 
fl <- list.files("~/ftc_atac/results/bwa/mergedLibrary/", full.names = TRUE)
fl <- fl[grep(".bam$", fl)]
# Let's quantify one file and see how long it takes...
Sys.time()
fc <- featureCounts(fl,
                    annot.ext = "~/ftc_atac/results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.saf",
                    fracOverlap = 0.2,
                    useMetaFeatures=TRUE,
                    isPairedEnd = TRUE,
                    autosort = FALSE,
                    # might want to increase this as well. 
                    nthreads = 6)
save(fc, file = "fc.RData")
Sys.time()