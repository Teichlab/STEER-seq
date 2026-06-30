## Needs to be run within the multiome environment

## conda activate multiome
## python 006_submit_Amulet.py $MANIFEST ($RUN_ID:optional)

import os
import sys
import pandas as pd


# CSV file with specific column names that will be querried by these scripts
## See example
## path2manifest = '/nfs/team205/is10/projects/thygenereg/analysis/001_preprocessing/thygenereg_by_sample_manifest - teichlab_new.csv'
path2manifest = sys.argv[1]

## Optionally specify 0-based list of indices to include (e.g., 0,1,2 will run for the first 3 samples)
if len(sys.argv) > 2:
    indices = [int(i) for i in sys.argv[2].split(',')]
else:
    indices = None
    

manifest = pd.read_csv(path2manifest)

if indices is None:
    indices = range(manifest.shape[0])

for i in indices:
    ## Make output directories
    os.makedirs((manifest['multiQC_dir'].iloc[i]+ '/amulet/'), exist_ok=True)

    ## Create single cell csv in required format
    barcodes = pd.read_csv((manifest['multiQC_dir'].iloc[i]+ '/multiQC_barcode_info.tsv'))
    barcodes_df = pd.DataFrame({'barcode':barcodes['barcode'], 'is__cell_barcode':1})
    barcodes_df.to_csv((manifest['multiQC_dir'].iloc[i]+ '/amulet/singlecell.csv'))

for i in range(manifest.shape[0]):
    if i in indices:
        ## Submit job
        os.system('/nfs/casm/team215mg/users/od5/scripts/jsub lsf -n amulet_' + manifest['sample_ID'].iloc[i] + ' -q normal -m 20g -c 8 -l logs "/nfs/casm/team215mg/users/od5/scripts/AMULET/AMULET.sh ' +  manifest['cellranger_path'].iloc[i] + '/atac_fragments.tsv.gz ' +  manifest['multiQC_dir'].iloc[i]+ '/amulet/singlecell.csv /nfs/casm/team215mg/users/od5/scripts/AMULET/human_autosomes.txt /nfs/casm/team215mg/users/od5/resources/hg38_repeatmasker/hg38_rmsk.bed ' + manifest['multiQC_dir'].iloc[i]+ '/amulet/ /nfs/casm/team215mg/users/od5/scripts/AMULET/" | bsub -G team215-grp')
