#!/bin/bash
set -eo pipefail

#map exactly like in N03, just park it into a PAF for downstream simplicity
~/minimap2/minimap2 -x map-hifi --secondary=no ~/transcriptomes+genemaps/GRCh38-2020-A/cdses-genenames.fa reads.fastq > mapped.paf
