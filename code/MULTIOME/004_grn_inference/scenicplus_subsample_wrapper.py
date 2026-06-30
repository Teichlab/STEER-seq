import dill
import pickle
import pandas as pd
import numpy as np
import argparse
from scenicplus.wrappers.run_scenicplus import run_scenicplus
import os

### This is a wrapper script to run scenicplus on a given initialised scenicplus object.
### is10 - teichlab, 01.2024

## Hardcoded annotations:
## Our annotation is based on Ensembl 98 so we can use the respective biomart host:
# biomart_host = "http://sep2019.archive.ensembl.org/" not available
biomart_host = "http://apr2020.archive.ensembl.org/"
## List of TF names
tf_file = '/nfs/casm/team215mg/users/od5/resources/utoronto_human_tfs_v_1.01.txt'
## Path to bigBed binary
path2bigBed = '/nfs/casm/team215mg/users/od5/scripts/binaries'
## Tmp dir for ray spilling
tmp_dir = '/lustre/scratch125/casm/teams/team215/r/'

## Parse arguments
parser = argparse.ArgumentParser(
                    prog='Scenicplus wrapper',
                    description='This is a wrapper script to run scenicplus on a given initialised scenicplus object.')

## Parse arguments from call
parser.add_argument('-w', '--work_dir', type=str, help='Directory - should contain an untrained or partially trained scenicplus object (naming scplus_obj.pkl) - all output files will be stored there as well.', required=True)
parser.add_argument('-c', '--ncpu', type=int, help='Number of threads to use', default=20, required=True)
parser.add_argument('-g', '--grouping_variable', type=str, help='Grouping variable (e.g., cell type)', required=True)
parser.add_argument('-s', '--seed', type=int, help='Random seed for subsampling cells', required=True)
parser.add_argument('-f', '--fraction', type=float, help='Fraction of cells to subsample', required=True)


args = parser.parse_args()

work_dir = args.work_dir
ncpu = args.ncpu
grouping_variable = args.grouping_variable
seed = args.seed
fraction = args.fraction

if fraction < 0 or fraction > 1:
    raise ValueError('fraction should be between 0 and 1')

## Load the full object
scplus_obj = pickle.load(open((work_dir + '/scplus_obj.pkl'), 'rb'))

## Make output directory
outdir = (work_dir + '/seed_' + str(seed) + '/')

## Make dir if it doesn't exist
os.makedirs(outdir, exist_ok=True)

## Subset
np.random.seed(seed)
cells = np.random.choice(scplus_obj.metadata_cell.index,
                         int(len(scplus_obj.metadata_cell.index) * fraction),
                         replace=False)
scplus_obj.subset(cells=cells)
pd.Series(cells).to_csv((outdir + '/cells.txt'), index=None, header=None)

## Run wrapper
try:
    run_scenicplus(
            scplus_obj = scplus_obj,
            variable = [grouping_variable],
            species = 'hsapiens',
            assembly = 'hg38',
            tf_file = tf_file,
            save_path = outdir,
            biomart_host = biomart_host,
            upstream = [1000, 150000],
            downstream = [1000, 150000],
            calculate_TF_eGRN_correlation = False,
            calculate_DEGs_DARs = False,
            export_to_loom_file = False,
            export_to_UCSC_file = False,
            path_bedToBigBed = path2bigBed,
            n_cpu = ncpu,
            _temp_dir = tmp_dir
    )
except Exception as e:
    #in case of failure, still save the object
    with open((outdir + "scplus_obj.pkl"), "wb") as f:
        dill.dump(scplus_obj, f)
    raise(e)

## Export unfiltered eregulons
scplus_obj.uns['eRegulon_metadata'].to_csv((outdir + 'scplus_eRegulon_metadata.csv'))

## Export unfiltered region to gene and TF2G links
scplus_obj.uns['region_to_gene'].to_csv((outdir + 'scplus_R2G.csv'))
scplus_obj.uns['TF2G_adj'].to_csv((outdir + 'scplus_TF2G.csv'))

## Filtering
from scenicplus.preprocessing.filtering import apply_std_filtering_to_eRegulons
apply_std_filtering_to_eRegulons(scplus_obj)

## Export filtered eregulons
scplus_obj.uns['eRegulon_metadata_filtered'].to_csv((outdir + 'scplus_eRegulon_metadata_filtered.csv'))

## Save object
with open((outdir + "scplus_obj.pkl"), "wb") as f:
    dill.dump(scplus_obj, f)