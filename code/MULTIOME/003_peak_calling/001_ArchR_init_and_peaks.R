## Should be run from within the archr conda environment
## conda activate archr
## /nfs/casm/team215mg/users/od5/scripts/jsub lsf -n archr -q normal -c 16 -m 250g "Rscript 001_ArchR_init_and_peaks.R" | bsub -G team215-grp

## Load required libraries
library(ArchR)
library(tidyverse)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)

## Initialised number of threads
addArchRThreads(threads = 30) 

## Genome and gene annotation
addArchRGenome("Hg38")

#### Parameters specific to the dataset - adjust appropriately 
## Path to arrow files
arrow_path <- '/lustre/scratch125/casm/teams/team215/users/od5/projects/celltherapy/ipsc/multiome/arrow_files/'

## Identify arrow files
ArrowFiles <- list.files(arrow_path, pattern='arrow$')
## Include only iPSC sample 
ArrowFiles <- ArrowFiles[grepl('Imm_NK15385635', ArrowFiles)==T]

cell_info_csv <- '/lustre/scratch125/casm/teams/team215/users/od5/projects/celltherapy/ipsc/multiome/gex_adata/tf_pilot/batch1/ipsc/leiden_v1/adata_multi_ipsc_perturbed_anno_v1_nuc_gex.csv'

## Specify path to macs2 (for peak calling)
pathToMacs2 <- '/software/conda/users/od5/archr_macs2/bin/macs2'

work_dir <- '/lustre/scratch125/casm/teams/team215/users/od5/projects/celltherapy/ipsc/multiome/archr/tf_pilot/batch1/ipsc/'
dir.create(work_dir, recursive=T)

## Set working directory (this is where the arrow files will be created)
setwd(work_dir)

## Initialise project
proj <- ArchRProject(
  ArrowFiles = paste0(arrow_path, ArrowFiles), 
  outputDirectory = "proj_immnk_tfpilot_ipsc",
  copyArrows = TRUE
)

## Add barcode name in same format as in annotated object to be able to transfer information
proj$sample_barcode <- paste0(proj$Sample, '-', gsub('.*#', '', proj$cellNames))

## Incorporate info about filtered dataset

cell_info <- read_csv(cell_info_csv)

archr_info <- getCellColData(proj) %>%
  as.data.frame() 

write_csv(archr_info, paste0(work_dir, 'archr_barcode_info.csv'))

archr_info <- archr_info %>%
  left_join(cell_info)

print(sum(proj$sample_barcode == archr_info$sample_barcode))

## Add some metadata to ArchR
cols2add <- c('sample_ID', 'sample_description', 'anno_v1_nuc_gex', 'perturbation')

lapply(cols2add, function(x) {
    proj <<- addCellColData(ArchRProj = proj, data = as.character(archr_info[, x]), name = x, cells = proj$cellNames)
})

## Save project
saveArchRProject(proj)

print(table(is.na(proj$anno_v1_nuc_gex)))

## Subset to cells with an annotation (i.e., those that passed qc)
proj <- subsetCells(proj,
                    cellNames = proj$cellNames[is.na(proj$anno_v1_nuc_gex)==FALSE])

## Adding group coverages
proj <- addGroupCoverages(
    ArchRProj = proj,
    groupBy = "perturbation", ## Group by cell type
    useLabels=FALSE, ## Don't use sample-aware pseudobulking as no replicate samples in experiment
    sampleLabels ="perturbation", ## And by perturbation
    minReplicates = 2, ## Minimum number of replicates (donors) required to include peak
    maxReplicates = 1000, ## We want to use all available replicates so we put this to an unreasonably high number
    minCells = 50, ## Minimum number of cells per replicate -> thus minimum number of cells per group will be 2 reps x 50 cells = 100 cells
    maxCells = 500, ## Maximum number of cells per replicate - Capping contribution to minimise bias of very abundant cell types
    sampleRatio = 0.8, ## Fraction of cells that can be sampled from a replicate to contribute to other replicates (without replacement)
    force=T
)

## And call peaks
proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "perturbation", 
    pathToMacs2 = pathToMacs2,
    reproducibility = "1", ## Minimum number of replicates (donors) a peak should be found in
    peaksPerCell = 5000, ## How many peaks do we allow per cell? -> i.e., a cell type with 50 cells can have max 50 X N peaks - we are lenient here and can filter downstream
    maxPeaks = 200000, ## Maximum number of peaks per cell type -> how many peaks can be identified per group? - we are lenient here and can filter downstream
    minCells = 5, ## The minimum allowable number of unique cells that was used to create the coverage files on which peaks are called (in practice, 80% of minCells in group coverages)
    excludeChr = c("chrY", "chrMT"), ## Chromosomes to exclude from peak calling
    method = "q", ## Filter by q-value
    cutOff = 0.01, ## Making sure we don't get very noisy peaks
    extendSummits = 250,  ## Extend summits by +/- 250 bp
    force = T
)

## Add peak matrix
proj <- addPeakMatrix(proj,
                      ceiling = 5,
                      binarize = F,
                      force = T)

## Save again
saveArchRProject(proj, outputDirectory = 'proj_immnk_tfpilot_ipsc')

## Export peak annotation

peaks <- getPeakSet(proj)
print(head(peaks))
print(length(peaks))

saveRDS(peaks, paste0(work_dir, "all_peaks_GRanges.rds"))

## And make bed file
peaks_bed <- data.frame(
    chr = seqnames(peaks),
    start = start(peaks),
    end = end(peaks),
    name = paste0(seqnames(peaks), ":", start(peaks), "-", end(peaks)) ## Peak name in UCSC format
)

write.table(peaks_bed, paste0(work_dir, "all_peaks.bed"), col.names=F, row.names=F, quote=F, sep='\t')

## Finally, getting the peak info (genomic context, peak calling scores etc)
peaks_info <- as.data.frame(peaks@elementMetadata)
peaks_info$name <- paste0(seqnames(peaks), ":", start(peaks), "-", end(peaks))

peaks_info <- left_join(peaks_info, peaks_bed)
write_tsv(peaks_info, paste0(work_dir, "all_peaks_info.txt"))

## Extracting peak matrix
## Summarized experiment
peak_se <- getMatrixFromProject(proj, "PeakMatrix")

## Actual matrix
peak_names <- paste0(seqnames(peak_se@rowRanges), ":", start(peak_se@rowRanges), "-", end(peak_se@rowRanges)) ## Peak names to match our bed files
peak_mat <- assay(peak_se)
row.names(peak_mat) <- peak_names

## Export peak matrix - we are using matrix format so that we can load in python as well
dir.create(paste0(work_dir, 'all_peaks_matrix'), recursive=T)
writeMM(peak_mat, paste0(work_dir, 'all_peaks_matrix/peak_mat.mtx'))
write(row.names(peak_mat), paste0(work_dir, 'all_peaks_matrix/peak_names.txt'))
write(colnames(peak_mat), paste0(work_dir, 'all_peaks_matrix/cell_names.txt'))

print('Done')
