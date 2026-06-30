## Needs to be run from within ArchR's R environment: 
## conda activate archr
## Rscript archr_makeArrows.R sample_id fragment_files output_p

## Script for making ArrowFile (for one sample)
## Requires as input the path to the fragment files, the sample ID and the path to the output file

## Currently hardcoded 
## TSS enrichment cutoff, minimum and maximum fragments per 

args = commandArgs(trailingOnly=TRUE)

sample_ID = args[1]
fragments = args[2]
output_p = args[3] 

## Load required libraries
library(ArchR)
library(tidyverse)
library(parallel)

## Initialised number of threads
addArchRThreads(threads = 8)

## Genome and gene annotation
addArchRGenome("Hg38")

## Set working directory (this is where the arrow files will be created)
setwd(output_p)

## Make arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = fragments,
  sampleNames = sample_ID,
  minTSS = 5, ## Minimum TSS enrichment (per cell)
  minFrags = 1e3, ## Minimum number of fragments (per cell)
  maxFrags = 1e5, ## Maximum number of fragments (per cell)
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  QCDir = "qc_dir",
  promoterRegion = c(2000,100),
  excludeChr = "chrY",
  force = T,
  subThreading = F, 
  cleanTmp = FALSE
)

## Done!