import os
import sys
import tables
import numpy as np
import scanpy as sc
import pandas as pd
import subprocess

## Running cellbender on a given sample 

#this is required so cellbender can be called from within the thing
os.environ['MKL_THREADING_LAYER'] = 'GNU'

cellranger_p=sys.argv[1]
cellbender_p=sys.argv[2]

## Directory should already have been made from within the submitting script

## First, read raw matrix and subset for GEX only
adata_raw = sc.read_10x_h5(cellranger_p + "/raw_feature_bc_matrix.h5")
adata_raw_gex = adata_raw[:,adata_raw.var.feature_types == 'Gene Expression']
adata_raw_gex.write_h5ad((cellbender_p + '/raw_gex_input.h5ad'))

## Read filtered matrix (from cellranger) and get the top 60% as very confident cells:
adata = sc.read_10x_h5(cellranger_p + "/filtered_feature_bc_matrix.h5")

## Run cellbender 
## Regarding priors, we want to be very confident that everything left of expected cells is cells and that everything right of the total droplets is empty - therefore, it makes sense to be quite conservative with these estimates
## We take the 60% of cells in the cellranger object as definite cells
## And we add 2 x nCells (or up to 15,000) barcodes to those of cellranger to reach definitely empty droplets
extra_barcodes = 2 * adata.obs.shape[0]
if extra_barcodes > 15000:
    extra_barcodes = 15000

def run_command(command):
    try:
        # Execute the command and capture stdout and stderr
        result = subprocess.run(command, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return result.stdout  # Normal output if command succeeds
    except subprocess.CalledProcessError as e:
        # This catches errors and returns the stderr where the error message is likely to be
        return e.stderr


## Construct base cellranger command
base_command = [
    "cellbender", "remove-background",
    "--input", f"{cellbender_p}/raw_gex_input.h5ad",
    "--output", f"{cellbender_p}/cellbender.h5",
    "--expected-cells", str(int(adata.obs.shape[0] * 0.6)),
    "--total-droplets-included", str(adata.obs.shape[0] + extra_barcodes),
    "--fpr", "0.01", "0.05", "0.1",
    "--cuda"
]

## Run the command (and catch any errors)
error_output = run_command(base_command)
if "AssertionError: There are no empty droplets with UMI counts over the lower cutoff" in error_output:
    print("Caught specific AssertionError, adjusting command")
    # Modify command if specific error caught
    base_command.extend(["--low-count-threshold", "1"])
    # Rerun the command with adjusted parameters
    final_output = run_command(base_command)
    print(final_output)
else:
    print(error_output)  # Print error message if other errors are caught
        
#fix the cellbender output files
#after this operation they should be importable via scanpy
for suf in ['_FPR_0.01_filtered.h5', '_FPR_0.01.h5', '_FPR_0.05_filtered.h5', '_FPR_0.05.h5', '_FPR_0.1_filtered.h5', '_FPR_0.1.h5']:
    h5_to_fix = cellbender_p+ '/cellbender' + suf
    with tables.open_file(h5_to_fix, "r+") as f:
        n = f.get_node("/matrix/features")
        n_genes = f.get_node("/matrix/shape")[0]
        if "genome" not in n:
            f.create_array(n, "genome", np.repeat("GRCh38", n_genes))