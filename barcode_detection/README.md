# Identifying transgene presence in sequencing data

Quantifying transgenes should be done based on the presence of the corresponding synthetic barcode. If you've created an additional library specifically amplifying the barcodes, then you can make use of Cellranger's feature reference functionality to get a count matrix.

If you're trying to extract the presence of transgene barcodes from regular gene expression libraries, there are two options:
- Explicitly search the sequence of the reads for the transgene barcodes. This can be accomplished with `parse_barcodes.py`, which depends on pysam and pandas. Provide the input BAM to `--bam`, provide a two-column CSV with the first column being barcode sequences and the second column being the corresponding transgene to `--seq`, and provide the output file name to `--csv`. The result is a CSV file with a count matrix for the identified transgenes in the identified cells. The script makes use of `CB`/`UB` tags in the BAM. It's recommended to pad the barcode sequences to 20bp with surrounding regions to avoid spurious matching.
- Create FASTA/GTF entries for the exogenous sequences, concatenate them with your existing reference's FASTA/GTF and use the extended reference for mapping. This approach requires careful construction of the exogenous sequences to avoid creating redundant regions - any read mapping to an exogenous sequence entry should clearly originate from that exact sequence. Practical considerations:
    - If the barcode is near (less than a read's length) from the transgene body, it's fine to include some of the transgene body in the exogenous sequence. However, avoid extending it beyond a read length, as then it becomes impossible to differentiate whether a read originates from the transgene or the native gene copy.
    - If the barcode is surrounded by long (more than a read's length) synthetic regions, provide the barcode surrounded by sufficiently short flanking synthetic sequences to force the mapped reads to include the barcode sequence.
    - The utility script `exog-fa-to-gtf.py` can be used to create GTF entries for exogenous sequence in `exog.fa`, populating requisite fields for Cellranger/Starsolo based on the name of the FASTA records. Here's an example FASTA record and the resulting GTF entry:
```
>EXOG_SPI1
AGAAGCCCgccaccatggggGGATCCTTACAGGCGTGCAAAATGGAAGGGTTTCCCCTCGTCCCCCCTCCATCAGAGGACCTGGTGCCCTATGACACGGATCTATACCAACGCCA
```
```
EXOG_SPI1       RefSeq  exon    1       115     .       +       .       transcript_id "EXOG_SPI1"; gene_id "EXOG_SPI1"; gene_name "EXOG_SPI1"
```

## Large-scale transgene screen

For the large-scale transgene screen, the information for which barcode corresponded to which transgene was not available a priori. The experimentally obtained constructs were sequenced via PacBio to get sufficiently long reads to guarantee the presence of both the barcode and the transgene body. The reads were then mapped to a set of human CDS sequences via the appropriate `minimap2` preset, obtaining alignment information in a PAF file. This PAF file is subsequently used by `find-barcodes-clean-paf.py` to take the start of the transgene body and the strandedness of the alignment, and mine the sequence upstream of it for the barcode via fuzzy matching flanking sequences. The identified barcode-transgene list is subsequently filtered to require the barcode to be between 10 and 14bp long, and for there to be more than 10 reads present. This cleans up a number of sequencing errors and lower evidence pairs.

Scripts for the analysis can be found in the subfolder.
