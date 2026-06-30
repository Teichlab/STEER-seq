## Should be run from within the archr conda environment
## conda activate archr
## /nfs/team205/is10/scripts/jsub lsf -n archr_qc -q long -c 16 -m 100g "Rscript 007_ArchR_TSSenr_FragPeriod.R" | bsub -G team205

suppressPackageStartupMessages({
library(ArchR)
library(tidyverse)
})
    
addArchRThreads(threads = 14) 

work_dir <- '/lustre/scratch126/cellgen/team205/is10/thygenereg/analysis/archr/'

## Load project
proj <- loadArchRProject(paste0(work_dir, 'proj_thymus/'))

## Subset to cells in anno_v3
cell_anno <- read_csv('/lustre/scratch126/cellgen/team205/is10/thygenereg/analysis/anno_v3/anno_v3_full.csv', guess_max =100000)

## Subset cells and add anno_v3
proj <- subsetCells(proj, cellNames = proj$cellNames[proj$sample_barcode %in% cell_anno$sample_barcode])

fragment_sizes <- plotFragmentSizes(ArchRProj = proj, returnDF = T)

head(fragment_sizes)
dim(fragment_sizes)

fragment_sizes <- as.data.frame(fragment_sizes) %>%
  dplyr::rename(Sample=group)

saveRDS(fragment_sizes, paste0(work_dir, "atac_QC_fragment_sizes_bySample.rds"))

p <-ggplot(fragment_sizes, aes(x=fragmentSize, y=fragmentPercent, color=as.factor(Sample))) +
  geom_line() +
  theme_classic() +
  xlab("Fragment size") +
  ylab("Density") +
  theme(legend.position = "bottom")+
  scale_color_manual(values = colorRampPalette(ArchRPalettes$stallion)(length(unique(fragment_sizes$Sample))), name="Sample")

pdf(paste0(work_dir,"atac_FragLength_distribution_byLib.pdf"), width = 6, height = 6, useDingbats = F); print(p);dev.off()

tssEnr <- plotTSSEnrichment(ArchRProj = proj, returnDF = T)

head(tssEnr)
dim(tssEnr)

tssEnr <- as.data.frame(tssEnr) %>%
  dplyr::rename(Sample=group) 

saveRDS(tssEnr, paste0(work_dir, "atac_QC_TSS_enrichment_bySample.rds"))

p <-ggplot(tssEnr, aes(x=x, y=smoothValue, color=as.factor(Sample))) +
  geom_line() +
  theme_classic() +
  xlab("Distance from TSS (bp)") +
  ylab("Normalised signal") +
  theme(legend.position = "bottom")+
  scale_color_manual(values = colorRampPalette(ArchRPalettes$stallion)(length(unique(tssEnr$Sample))), name="Sample")

pdf(paste0(work_dir, "atac_TSS_signal_distribution_byLib.pdf"), width = 6, height = 6, useDingbats = F); print(p);dev.off()

print('Done!')