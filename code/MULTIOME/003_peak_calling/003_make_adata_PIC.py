## conda activate scvi_24
## /nfs/casm/team215mg/users/od5/scripts/jsub lsf -n make_adata_pic -q normal -c 1 -m 500g "python 003_make_adata_PIC.py" | bsub -G team215-grp 

import os
import sys
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata as ad
import scvi
import session_info

from scipy.io import mmread
from scipy.sparse import csr_matrix

sys.path.append('/nfs/casm/team215mg/users/od5/scripts/')
from single_cell_utils import *
from scjp_colors import *

## Set up work_dir
work_dir = '/lustre/scratch125/casm/teams/team215/users/od5/projects/celltherapy/ipsc/multiome/archr/tf_pilot/batch1/ipsc/pic_mat/'

## Make PIC object
## Read matrices and make adata
input_dir = work_dir  + 'by_sample/'
samples = ['multi_Imm_NK15385635']

adata_list = []
for s in samples:
    print(s)
    mat = csr_matrix(mmread((input_dir + s + '/peak_mat.mtx')))
    
    cell_names = pd.read_csv((input_dir + s + '/cell_names.txt'), header=None)[0]
      
    peak_names = pd.read_csv((input_dir + s + '/peak_names.txt'), header=None)[0]

    adata = ad.AnnData(X=mat.transpose())
    adata.obs_names = cell_names
    adata.var_names = peak_names

    adata_list.append(adata)

## Put everything together
adata_full = ad.concat(adata_list,join='outer')

## Export adata obs
adata_full.write_h5ad('/lustre/scratch125/casm/teams/team215/users/od5/projects/celltherapy/ipsc/multiome/archr/tf_pilot/batch1/ipsc/pic_mat/pic_adata.h5ad')
