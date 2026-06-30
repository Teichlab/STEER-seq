## Call this from the multiome environment

## conda activate multiome
## python 003_submit_multiQC.py $MANIFEST ($indices:optional)

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
    
## Path to the script multiQC_single_sample.ipynb
script_dir = '/nfs/casm/team215mg/users/od5/projects/celltherapy/ipsc/multiome/tf_pilot/001_preprocessing/'

manifest = pd.read_csv(path2manifest)

if indices is None:
    indices = range(manifest.shape[0])

for i in range(manifest.shape[0]):
    if i in indices:
        ## Make output directory
        os.makedirs(manifest['multiQC_dir'].iloc[i], exist_ok=True)
        ## Copy script in location
        os.system('cp ' + script_dir + "multiQC_single_sample.ipynb " + manifest['multiQC_dir'].iloc[i])
        ## Change directory to output directory
        os.chdir(manifest['multiQC_dir'].iloc[i])
        ## Submit job
        os.system('/nfs/casm/team215mg/users/od5/scripts/jsub lsf -n ' + manifest['sample_ID'].iloc[i] + ' -q normal -m 50g -l logs "papermill multiQC_single_sample.ipynb multiQC_single_sample_' + manifest['sample_ID'].iloc[i] + '.ipynb -p sample ' + manifest['sample_ID'].iloc[i] + ' -p input_dir ' + manifest['cellranger_path'].iloc[i] + ' -p output_dir ' + manifest['multiQC_dir'].iloc[i] + ' -p cell_calling_mode ' + manifest['cell_calling_mode'].iloc[i] + ' -p topN_cells ' + str(manifest['topN_cells'].iloc[i]) + ' -p celltypist_mod ' + manifest['celltypist_mod'].iloc[i] + ' -p gex_median_umi_cutoff ' + manifest['gex_median_umi_cutoff'].iloc[i] + ' -p emptyDropsMultiOut ' + manifest['emptyDropsMulti_dir'].iloc[i] + ' -p markerGenes ' + manifest['markerGenes'].iloc[i] + '" | bsub -G team215-grp')
