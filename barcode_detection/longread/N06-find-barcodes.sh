#!/bin/bash
set -eo pipefail

python3 find_barcodes_clean_paf.py reads.fastq mapped_cleaned.paf barcodes.csv
