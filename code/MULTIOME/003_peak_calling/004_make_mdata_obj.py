## conda activate scvi_24
## /nfs/casm/team215mg/users/od5/scripts/jsub lsf -n make_mdata -q normal -c 1 -m 200g "python 004_make_mdata_obj.py" | bsub -G team215-grp 

import os
import sys
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata as ad
import session_info
import muon as mu
from muon import MuData

sys.path.append('/nfs/casm/team215mg/users/od5/scripts/')
from single_cell_utils import *
from scjp_colors import *

work_dir = '/lustre/scratch125/casm/teams/team215/users/od5/projects/celltherapy/ipsc/multiome/'

## Make mdata object

## Read RNA and ATAC
rna = sc.read_h5ad(work_dir + 'gex_adata/tf_pilot/batch1/ipsc/leiden_v1/adata_leiden_v3_ipsc_perturbed.h5ad')
atac = sc.read_h5ad(work_dir + 'archr/tf_pilot/batch1/ipsc/pic_mat/pic_adata.h5ad')

## Subset to common cells
common_cells = rna.obs_names[rna.obs_names.isin(atac.obs_names)]

## Cells to be excluded:
print(rna.obs['sample_ID'][rna.obs_names.isin(common_cells)==False].value_counts().head(10))
print(rna.obs['perturbation'][rna.obs_names.isin(common_cells)==False].value_counts().head(10))

rna = rna[common_cells]
atac = atac[common_cells]

## Add rna obs to ATAC for convenience
atac.obs = rna.obs.copy()

mdata = MuData({"rna": rna, "atac": atac})

## Export for further use:
mdata.write_h5mu((work_dir + 'mdata/tf_pilot/batch1/ipsc/mdata_full.h5mu'))
