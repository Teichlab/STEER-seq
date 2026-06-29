#!/bin/bash
set -eo pipefail

#we just have one read file this time
echo "mapped,total" > summary.csv

~/minimap2/minimap2 -ax map-hifi --secondary=no ~/transcriptomes+genemaps/GRCh38-2020-A/cdses-genenames.fa reads.fastq > mapped.sam

#the sam has three primary types of records
#a flag of 2048 means a "supplementary" (chimeric) alignment
#in the even of barcode PAF, those reads would have multiple entries
#and can potentially be discarded on the basis of duplicate indices
RCOUNT=`samtools view -F 2048 mapped.sam | wc -l`
#the second kind of entry has a flag of 4 for unmapped
#ignore those too and you get the third kind, the regular mapped
MCOUNT=`samtools view -F 4 -F 2048 mapped.sam | wc -l`
#all the info we need for the report
echo "${MCOUNT},${RCOUNT}" >> summary.csv
