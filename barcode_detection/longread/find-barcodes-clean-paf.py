from rapidfuzz.fuzz import ratio, partial_ratio_alignment
import pandas as pd
import time
import sys

#TODO: make this nicer with argparse but can't be bothered right now
read_path = sys.argv[1]
paf_path = sys.argv[2]
output_path = sys.argv[3]
pre_barcode = "CGCGCCGAAGCGACC"
post_barcode = "CTTGTCTCGCTCCGG"
score_thresh = 90

def revcomp(seq):
    return seq.translate(str.maketrans('ATGC','TACG'))[::-1]


#stash outputs in a dictionary as we go
hits = {}

#cleaned up paf via clean_paf.py
#will read from it iteratively as we go
fid2 = open(paf_path, "r")
#skip header
line = fid2.readline()
#nab first read with info
line = fid2.readline()
#0 - start, 1 - end, 2 - orientation, 3 - target
paf_name = line.split(",")[0]
paf_single = line.rstrip().split(",")[1:]
#the first two are integers
paf_single[0] = int(paf_single[0])
paf_single[1] = int(paf_single[1])

#may as well time it because why not
t1 = time.time()
with open(read_path, "r") as fid:
    lines = []
    for line in fid:
        lines.append(line.rstrip())
        if len(lines) == 4:
            #split on spacebar and lose the leading @
            name = lines[0].split(" ")[0][1:]
            if name == paf_name:
                sequence = lines[1]
                #find the part that can have the barcode
                if paf_single[2] == "+":
                    #everything up to start is fair game
                    sequence = sequence[:paf_single[0]]
                else:
                    #everything after end is fair game, but needs a revcomp
                    sequence = revcomp(sequence[paf_single[1]:])
                #find barcode flanking sequence hits
                pra_pre = partial_ratio_alignment(sequence, pre_barcode)
                #require that both alignments be decent quality
                #and that the two sequences are identified in order
                #check stuff on the fly for run time purposes
                #nested if bullshit time!
                if pra_pre.score >= score_thresh:
                    #do a search relative to the end of the previous one
                    pra_post = partial_ratio_alignment(sequence[pra_pre.src_end:], post_barcode)
                    if (pra_post.score > score_thresh):
                        #oh nice a barcode! don't forget the second ratio starts where the first ends
                        barcode = sequence[pra_pre.src_end:(pra_pre.src_end+pra_post.src_start)]
                        #stash the barcode and target, and the two scores for science
                        #maybe the losses will be low enough to not matter?
                        hits[name] = [barcode, paf_single[3], pra_pre.score, pra_post.score]
                #need to extract new expected read from paf
                line = fid2.readline()
                #end of file will manifest as this being an empty string
                if line == "":
                    break
                paf_name = line.split(",")[0]
                paf_single = line.rstrip().split(",")[1:]
                #the first two are integers
                paf_single[0] = int(paf_single[0])
                paf_single[1] = int(paf_single[1])
            #blank holder for next upcoming record
            lines = []

t2 = time.time()
print(t2-t1)

#turn to df
hitsdf = pd.DataFrame.from_dict(hits, orient="index", columns=["barcode", "target", "pre_score", "post_score"])
hitsdf.to_csv(output_path)
