## Submit from within the multiome conda environment
#conda activate multiome
#/nfs/casm/team215mg/users/od5/scripts/jsub lsf -n pycistgt -q long -c 1 -m 200g "python 004_pycistarget_topics_and_DARs.py" | bsub -G team215-grp

import os
import pickle
import dill
import pandas as pd

import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
from scenicplus.wrappers.run_pycistarget import run_pycistarget

## Set up work directory
work_dir = '/lustre/scratch125/casm/teams/team215/users/od5/projects/celltherapy/ipsc/multiome/scenicplus/tf_pilot/batch1/ipsc/'

## Create directory to store the motifs:
os.makedirs(os.path.join(work_dir, 'motifs'), exist_ok=True)

## Tmp dir for ray (this needs to be as short as possible, as there is an overall limit of the number of characters accepted by ray - max is 107):
tmp_dir = '/lustre/scratch125/casm/teams/team215/r/'
os.makedirs(tmp_dir, exist_ok=True)

## Load binarised topic regions and DARs from previous step
region_bin_topics_otsu = pickle.load(open((work_dir  + '/cistopic/scATAC/region_bin_topics_otsu.pkl'), 'rb'))
region_bin_topics_top3k = pickle.load(open((work_dir  + '/cistopic/scATAC/region_bin_topics_top3K.pkl'), 'rb'))

dars_archr_df = pd.read_csv((work_dir  + '/DARs_archr/archr_DARs_perturbation_df.txt'), sep='\t')
dars_archr_df['Cell_type'] = dars_archr_df['Cell_type'].replace('/', '_', regex=True)

## The topics are already in the correct format - for the DARs, converting to dictionaries
dars_archr = {}
for ct in dars_archr_df['Cell_type'].unique():
    dars_archr[ct] = dars_archr_df['peak'][dars_archr_df['Cell_type']==ct]

## Convert to pyranges
## Note that we have a different structure in the topics/dars so different code is needed:
region_sets = {}
region_sets['topics_otsu'] = {}
region_sets['topics_top_3'] = {}
region_sets['DARs_archr'] = {}

for topic in region_bin_topics_otsu.keys():
    regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for topic in region_bin_topics_top3k.keys():
    regions = region_bin_topics_top3k[topic].index[region_bin_topics_top3k[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_top_3'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for DAR in dars_archr.keys():
    regions = dars_archr[DAR][dars_archr[DAR].str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['DARs_archr'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))

## Export this
pickle.dump(region_sets,
            open((work_dir + 'region_sets.pkl'), 'wb'))

## Load cistarget databases (for now using precomputed databases, we need to eventually shift to our own, probably with new peak calls)
db_fpath = '/lustre/scratch125/casm/teams/team215/users/od5/projects/celltherapy/ipsc/multiome/archr/tf_pilot/batch1/ipsc/cistarget_db/'
motif_annot_fpath = "/nfs/casm/team215mg/users/od5/resources/aerts_motifs/v10nr_clust_public/snapshots/"

rankings_db = os.path.join(db_fpath, 'cistarget_db.regions_vs_motifs.rankings.feather')
scores_db =  os.path.join(db_fpath, 'cistarget_db.regions_vs_motifs.scores.feather')
motif_annotation = os.path.join(motif_annot_fpath, 'motifs-v10-nr.hgnc-m0.00001-o0.0.tbl')

## Run pycistarget wrapper
run_pycistarget(
    region_sets = region_sets,
    species = 'homo_sapiens',
    save_path = os.path.join(work_dir, 'motifs'),
    ctx_db_path = rankings_db,
    dem_db_path = scores_db,
    path_to_motif_annotations = motif_annotation,
    run_without_promoters = True,
    n_cpu = 1,
    _temp_dir = tmp_dir,
    annotation_version = 'v10nr_clust',
    )

## Done
print("Done!")