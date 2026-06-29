import pandas as pd

paf_path="mapped.paf"

#PAF format: https://github.com/lh3/miniasm/blob/master/PAF.md
#we want the read name, the start/end/orientation on the query, and target name
paf = pd.read_table(paf_path, index_col=0, header=None, usecols=[0,2,3,4,5], names=["", "start", "end", "orientation", "target"])
#this identifies all duplicates, negate to extract all the single mapped ones
paf = paf.loc[~paf.index.duplicated(keep=False), :]
paf.to_csv("mapped_cleaned.paf")
