A collection of scripts and notebooks capturing the process of reconstructing barcode-transgene relationships from PacBio reads

1. Getting the reads into `reads.fastq` (skipped, pulling from internal iRODS)
2. Preparing a CDS FASTA file from 10X's GRCh38-2020-A reference
3. Mapping the data with minimap2 to get mapping rates
4. Mapping the data with minimap2 to get the alignments in a PAF file
5. Cleaning up the PAF file, removing duplicates and keeping a subset of columns
6. Invoking `find-barcodes-clean-paf.py` to identify barcodes by fuzzy matching flanking sequences, and reporting barcode-transgene combinations on a per-read basis
7. Cleaning up intermediate files (skipped)
8. Exploratory analysis of the barcode-transgene data, flagging the issue of occasional sequencing errors
9. Applying a 10-read cutoff to the data, handling the sequencing errors, and generating a Cellranger feature reference file
