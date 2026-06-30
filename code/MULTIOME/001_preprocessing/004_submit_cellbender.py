## Call this from the cellbender environment

## conda activate cellbender
## python 004_submit_cellbender.py $MANIFEST ($indices:optional)

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

## Path to the script cellbender_wrapper.py
script_dir = '/nfs/casm/team215mg/users/od5/projects/celltherapy/ipsc/multiome/tf_pilot/001_preprocessing/'

manifest = pd.read_csv(path2manifest)

if indices is None:
    indices = range(manifest.shape[0])

for i in range(manifest.shape[0]):
    if i in indices:
        ## Make output directory
        os.makedirs(manifest['cellbender_dir'].iloc[i], exist_ok=True)
        ## Submit job
        os.system(('bsub -J cellbender_' + manifest['sample_ID'].iloc[i] + ' -G team215-grp -q gpu-normal -n1 -M15000 -R"select[mem>15000] rusage[mem=15000] span[hosts=1]" -gpu"mode=shared:j_exclusive=no:gmem=10000:gmodel=NVIDIAA100_SXM4_80GB:num=1" -o logs/cellbender_'+ manifest['sample_ID'].iloc[i] +'.out -e logs/cellbender_' + manifest['sample_ID'].iloc[i] +'.err "python cellbender_wrapper.py ' + manifest['cellranger_path'].iloc[i] + ' ' + manifest['cellbender_dir'].iloc[i] +  '"'))
        
        #os.system('/software/team205/bin/jsub lsf -n cellbender_' + manifest['sample_ID'].iloc[i] + ' -q gpu-normal -c 1 -m 15g "python cellbender_wrapper.py ' + manifest['cellranger_path'].iloc[i] + ' ' + manifest['cellbender_dir'].iloc[i] +  '" | bsub -G teichlab') -> the problem with this sub is that we can't specify cuda memory (will just use what is available) and it will occasionally run out of memory
