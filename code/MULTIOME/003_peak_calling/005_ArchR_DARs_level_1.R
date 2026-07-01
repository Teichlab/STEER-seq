## Should be run from within the archr conda environment
## conda activate archr
## /nfs/casm/team215mg/users/od5/scripts/jsub lsf -n archr_dars_1 -q normal -c 16 -m 100g "Rscript 005a_ArchR_DARs_level_1.R" | bsub -G team215-grp

suppressPackageStartupMessages({
library(ArchR)
library(tidyverse)
})
    
addArchRThreads(threads = 30) 

work_dir <- '/lustre/scratch125/casm/teams/team215/users/od5/projects/celltherapy/ipsc/multiome/archr/tf_pilot/batch1/ipsc/'

## Load project
proj <- loadArchRProject(paste0(work_dir, 'proj_immnk_tfpilot_ipsc/'))

## Subset to cells in anno_v3
cell_anno <- read_csv('/lustre/scratch125/casm/teams/team215/users/od5/projects/celltherapy/ipsc/multiome/gex_adata/tf_pilot/batch1/ipsc/leiden_v1/adata_multi_ipsc_perturbed_anno_v1_nuc_gex.csv', guess_max =100000)

## Subset cells and add perturbation_cluster
proj <- subsetCells(proj, cellNames = proj$cellNames[proj$sample_barcode %in% cell_anno$sample_barcode])

archr_info <- getCellColData(proj) %>%
  as.data.frame() %>%
  left_join(select(cell_anno, sample_barcode, anno_v1_nuc_gex, perturbation))
print(nrow(archr_info))

cols2add <- c('perturbation')

lapply(cols2add, function(x) {
    proj <<- addCellColData(ArchRProj = proj, data = as.character(archr_info[, x]), name = x, cells = proj$cellNames, force = TRUE)
})

## Calling marker peaks based on the peak-groups
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "perturbation",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "ttest"
)

saveRDS(markersPeaks, paste0(work_dir, '/archr_DARs_perturbation_obj.rds'))

## Extract markers
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.2 & Log2FC >= 0")
markerList

markerPeaks_df <- Reduce(bind_rows,lapply(1:length(markerList), function(i) {
  ## Cell type or state
  ct <- names(markerList)[i]
  ## Object of marker peaks
  gr <- markerList[[i]]
  if (nrow(gr) == 0) {
    return(NULL)
  } else {
     ## Getting peak names and stats
    peaks <- data.frame(chr=as.character(gr$seqnames), 
                        start=gr$start,
                        end=gr$end,
                        Log2FC=gr$Log2FC,
                        FDR=gr$FDR,
                        MeanDiff=gr$MeanDiff,
                        Cell_type=ct,
                        stringsAsFactors = F) %>%
      mutate(peak=paste0(chr, ':',start, '-',end))
    return(peaks)
  }
}))

## Export
write_tsv(markerPeaks_df, paste0(work_dir, 'archr_DARs_perturbation_df.txt'))

## Plot heatmap
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.2 & Log2FC >= 0", 
  transpose = FALSE,
  binaryClusterRows = T, 
  clusterCols = T,
  nLabel = 3,
  nPrint = 3
)

pdf(paste0(work_dir, 'archr_DARs_perturbation_heatmap.pdf'), width =10, height = 8); draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot");dev.off()
