## Activate conda environment before running the script

## Example usage
#conda activate create_cistarget_databases
#/nfs/casm/team215mg/users/od5/scripts/jsub lsf -n cistarget_db -q long -c 16 -m 200g "bash 008_create_cistarget_db.sh" | bsub -G team215-grp

## Location of all peaks to be used as background for cistarget analyses
peaks_fasta='/lustre/scratch125/casm/teams/team215/users/od5/projects/celltherapy/ipsc/multiome/archr/tf_pilot/batch1/ipsc/all_peaks_plus_encode.fa'
## Directory of cistarget database to be created
out_dir='/lustre/scratch125/casm/teams/team215/users/od5/projects/celltherapy/ipsc/multiome/archr/tf_pilot/batch1/ipsc/cistarget_db'

mkdir -p $out_dir

## Set ${create_cistarget_databases_dir} variable to path where the repo was cloned to.
create_cistarget_databases_dir='/nfs/casm/team215mg/users/od5/scripts/create_cisTarget_databases'
motif_collection='/nfs/casm/team215mg/users/od5/resources/aerts_motifs/v10nr_clust_public'

## Create database
python ${create_cistarget_databases_dir}/create_cistarget_motif_databases.py \
-f $peaks_fasta \
-M ${motif_collection}/singletons/ \
-m ${motif_collection}/motifs.txt \
-o $out_dir \
-s 42 \
-t 14

## Create rankings
${create_cistarget_databases_dir}/convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbs.py \
-i $out_dir.motifs_vs_regions.scores.feather \
-s 42
