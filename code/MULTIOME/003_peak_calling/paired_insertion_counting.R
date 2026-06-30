## Needs to be run from within ArchR's R environment: 
## conda activate archr
## /nfs/team215/od5/scripts/jsub lsf -n pic_$sample_ID -q normal -c 4 -m 40g "Rscript paired_insertion_counting.R $sample_ID" | bsub -Gteam215-grp 

## Script for counting paired insertions (for one sample)
## Requires as input the path to the fragment files, the sample ID and the path to the output file

## Currently hardcoded 
## Path to peak set

args = commandArgs(trailingOnly=TRUE)

sample_ID = args[1]

print(sample_ID)

## Load libraries
library('data.table')
library('tidyverse')
library('GenomicRanges')
library('Matrix')
library('PICsnATAC')

## Output directory
output_p <- paste0('/lustre/scratch125/casm/teams/team215/users/od5/projects/celltherapy/ipsc/multiome/archr/tf_pilot/batch1/ipsc/pic_mat/by_sample/', sample_ID , '/')
if (dir.exists(output_p) == F) {
    dir.create(output_p, recursive=T)
}

## Read cell info
cell_info <-read_csv('/lustre/scratch125/casm/teams/team215/users/od5/projects/celltherapy/ipsc/multiome/gex_adata/tf_pilot/batch1/ipsc/leiden_v1/adata_multi_ipsc_perturbed_anno_v1_nuc_gex.csv')
# cell_info <- cell_info %>% filter(anno_v1_nuc_gex == 'iPSC')

## Read sample metadata 
metad <- read_csv('/nfs/casm/team215mg/users/od5/projects/celltherapy/ipsc/multiome/tf_pilot/001_preprocessing/immnk_by_sample_manifest_updated - immnk_20251028.csv')

## Combine and filter for samples in dataset
combined_metad <- metad %>%
filter(sample_ID %in% c("multi_Imm_NK15385635"))

## Extract raw barcodes and path to fragment files
frags_path <- paste0(combined_metad$cellranger_path[combined_metad$sample_ID==sample_ID], '/atac_fragments.tsv.gz')
cells <- cell_info$barcode[cell_info$sample_ID==sample_ID]

print(length(cells))
print(sum(is.na(cells)))

## Read peak set
peak_set <- readRDS('/lustre/scratch125/casm/teams/team215/users/od5/projects/celltherapy/ipsc/multiome/archr/tf_pilot/batch1/ipsc/all_peaks_GRanges.rds')

## Fix peak names
names(peak_set) <- paste0(peak_set@seqnames, ':', start(peak_set), '-', end(peak_set))

## Do the actual counting
pic_mat <- PIC_counting(
  cells = cells,
  fragment_tsv_gz_file_location = frags_path,
  peak_sets = peak_set
)

## Modify barcodes to "sample_barcode" format
colnames(pic_mat) <- paste0(sample_ID, '-', colnames(pic_mat))

## And export
writeMM(pic_mat, paste0(output_p, 'peak_mat.mtx'))
write(row.names(pic_mat), paste0(output_p, '/peak_names.txt'))
write(colnames(pic_mat), paste0(output_p, '/cell_names.txt'))

print('Done!')
