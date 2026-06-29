import numpy as np

#helper script to turn a transgene fasta into its corresponding gtf
#with one exon entry per sequence in the fasta

#load fasta into a single list of lines
with open("exog.fa", "r") as fid:
        exog = [i.rstrip() for i in fid.readlines()]

#turn into two-column numpy array with gene name in [:,0] and sequence in [:,1]
#the default reshape order works for this, no need to fiddle with anything
exog = np.reshape(exog, (int(len(exog)/2),2))

#build a GTF from the info we have
with open("exog.gtf", "w") as fid:
        for i in np.arange(exog.shape[0]):
                #strip out the ">" at the start
                gene = exog[i,0][1:]
                #we just care about the length, and in string form too
                length = str(len(exog[i,1]))
                #the last field of the GTF is a bunch of space delimited stuff
                #sneak in a newline at the end
                lastfield=" ".join(["transcript_id",'"'+gene+'";',"gene_id",'"'+gene+'";',"gene_name",'"'+gene+'"\n'])
                #the rest of the GTF line is tab joined
                fid.write("\t".join([gene,"RefSeq","exon","1",length,'.','+','.',lastfield]))
