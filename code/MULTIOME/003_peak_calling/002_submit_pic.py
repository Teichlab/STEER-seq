## Needs to be run from within the archr conda environment:
## conda activate archr
## python 002_submit_pic.py

import os
import pandas as pd

cell_info = pd.read_csv('/lustre/scratch125/casm/teams/team215/users/od5/projects/celltherapy/ipsc/multiome/gex_adata/tf_pilot/batch1/ipsc/leiden_v1/adata_multi_ipsc_perturbed_anno_v1_nuc_gex.csv', index_col=0)

for sample_ID in cell_info['sample_ID'].unique():
    os.system('/nfs/casm/team215mg/users/od5/scripts/jsub lsf -n pic_' + sample_ID + ' -q normal -c 4 -m 40g "Rscript paired_insertion_counting.R ' + sample_ID + '" | bsub -G team215-grp')
