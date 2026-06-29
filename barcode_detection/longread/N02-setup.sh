#downloaded minimap2 from github and compiled it
#prepared CDS-only file from gffread - as per Ivan, the synthetic sequences are CDSes
gffread -x ~/transcriptomes+genemaps/GRCh38-2020-A/cdses.fa -g ~/cellranger/GRCh38-2020-A/fasta/genome.fa ~/cellranger/GRCh38-2020-A/genes/genes.gtf
#this has a header of >ENST gene=GENENAME, turn it to >ENST_GENENAME
sed "s/ gene=/_/" ~/transcriptomes+genemaps/GRCh38-2020-A/cdses.fa > ~/transcriptomes+genemaps/GRCh38-2020-A/cdses-genenames.fa
