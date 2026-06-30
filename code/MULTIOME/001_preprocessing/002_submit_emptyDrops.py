## Call this from the emptyDrops environment

## conda activate eDm_env
## python 002_submit_emptyDrops.py $MANIFEST ($indices:optional)

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

for i in range(manifest.shape[0]):
    if i in indices:
        ## Submit emptyDrops job
        os.makedirs(manifest['emptyDropsMulti_dir'].iloc[i], exist_ok = True)
        os.system('/nfs/casm/team215mg/users/od5/scripts/jsub lsf -n eDm_' + manifest['sample_ID'].iloc[i] + ' -q normal -m 30g -l logs "Rscript emptyDrops.R ' + manifest['cellranger_path'].iloc[i] + '/raw_feature_bc_matrix ' +  manifest['emptyDropsMulti_dir'].iloc[i] + '" | bsub -G team215-grp')
