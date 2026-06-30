## Needs to be run within the archr environment

## conda activate archr
## python 005_submit_makeArrows.py $MANIFEST ($RUN_ID:optional)

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
    
## Path to the script archr_makeArrows.R
script_dir = '/nfs/casm/team215mg/users/od5/projects/celltherapy/ipsc/multiome/tf_pilot/001_preprocessing/'

manifest = pd.read_csv(path2manifest)

## Path to arrow files:
arrow_p  = '/lustre/scratch125/casm/teams/team215/users/od5/projects/celltherapy/ipsc/multiome/arrow_files/'
os.makedirs(arrow_p, exist_ok=True)

if indices is None:
    indices = range(manifest.shape[0])

for i in range(manifest.shape[0]):
    if i in indices:
        ## Submit job
        os.system('/nfs/casm/team215mg/users/od5/scripts/jsub lsf -n arrows_' + manifest['sample_ID'].iloc[i] + ' -q normal -m 40g -c 8 -l logs "Rscript archr_makeArrows.R ' + manifest['sample_ID'].iloc[i] + ' ' +  manifest['cellranger_path'].iloc[i] + '/atac_fragments.tsv.gz ' + arrow_p + '" | bsub -G team215-grp')
