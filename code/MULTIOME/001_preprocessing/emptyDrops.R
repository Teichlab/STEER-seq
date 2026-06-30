#!/usr/bin/env Rscript

## Script for running EmptyDrops
## Requires as input the path to the raw_feature_bc_matrix, as given by cellranger-arc

args = commandArgs(trailingOnly=TRUE)

input_p = args[1]
output_p = args[2] 


if (!require("BiocManager", quietly = TRUE))
	    install.packages("BiocManager", repos = "http://cran.us.r-project.org")

BiocManager::install("DelayedArray", force = TRUE)
BiocManager::install("S4Arrays", force = TRUE)

## Load libraries
library(Seurat)
library(EmptyDropsMultiome)

## Read matrix into sce object
sce <- Read10X(input_p)

## Split modalities into RNA and ATAC
rna <- sce[["Gene Expression"]]
atac <- sce[["Peaks"]]

## The main function - running emptydrops on the two matrices
pdf(paste0(output_p, "/plots.pdf"), width = 6, height=4, useDingbats = F)
eD.out <- emptydrops_multiome(count_matrix_rna=rna, count_matrix_atac=atac)
dev.off()

write.table(eD.out, paste0(output_p, "/emptyDropsMulti.csv"), col.names=T, row.names=F, sep=',', quote=F)
