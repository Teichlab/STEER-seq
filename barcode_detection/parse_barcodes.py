import argparse
import pandas as pd
import pysam
import os

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam', dest='bam', required=True, help='Mapped BAM file with transgene barcodes present; index (.bam.bai) frequired')
    parser.add_argument('--seq', dest='seq', required=True, help='CSV file with transgene barcodes to identify; first row header, first column barcode sequences, second column assigned barcode names')
    parser.add_argument('--tail', dest='tail', default=None, help='Sequence to append at the end of each specified barcode')
    parser.add_argument('--csv', dest='csv', default='barcode_counts.csv', help='Output CSV file to store the transgene barcode counts in')
    args = parser.parse_args()
    return args

def has_barcode(read, seqdict, rcseqdict):
    #this is atac so this can show up anywhere in the read
    for bc in seqdict:
        if bc in read.seq:
            return seqdict[bc]
    for rc in rcseqdict:
        if rc in read.seq:
            return rcseqdict[rc]
    return ''

def parse_barcodes(bam, seq):
    '''
    Create a data frame of cell barcode by multiplexing barcode counts.
    
    Input:
     * bam - path to an indexed BAM file with the multiplexing barcodes present at the 
     start of R2; index (.bam.bai) required
     * seq - pandas data frame with multiplexing barcode sequences as the index, and the 
     corresponding names in the first column
    '''
    #turn the df to a dictionary for ease of access, use the first column as the names
    seqdict = seq.to_dict()[seq.columns[0]]
    #create entries for reverse complements which will be at the end
    rcseqdict = {}
    for bc in seqdict:
        rc = bc.translate(str.maketrans('ATGC','TACG'))[::-1]
        rcseqdict[rc] = seqdict[bc]
    #we can pop open the BAM now
    bam_file = pysam.AlignmentFile(bam, "rb")
    #this will be a nested dictionary of lists - sequencing barcode, then cell barcode, then list of UMIs
    parsed = {}
    #track how many reads we hit that have the barcode, and matching cell/UMI information
    readcount = 0
    startswith = 0
    isproper = 0
    #don't skip unmapped - that's where most of the info resides!
    for read in bam_file.fetch(until_eof=True):
        readcount = readcount+1
        #parse the read - do we have a hit of one of the barcodes?
        bc = has_barcode(read, seqdict, rcseqdict)
        if bc != '':
            #alright, we found a read that has the barcode. but can we use it?
            startswith = startswith + 1
            if read.has_tag("CB") and read.has_tag("UB"):
                #we have a hit! time to write it
                cb = read.get_tag("CB")
                umi = read.get_tag("UB")
                #create required bits of nested dictionary 
                if bc not in parsed:
                    parsed[bc] = {}
                if cb not in parsed[bc]:
                    parsed[bc][cb] = []
                #stash UMI
                parsed[bc][cb].append(umi)
                isproper = isproper+1
    print("Saw "+str(readcount)+" reads")
    print("Found "+str(startswith)+" reads total with barcodes present")
    print("Of those, "+str(isproper)+" have cell information")
    #now we can replace the list of UMIs with the number of unique UMIs total
    for bc in parsed:
        for cb in parsed[bc]:
            parsed[bc][cb] = len(set(parsed[bc][cb]))
    #and now we can easily turn this to a data frame of counts and return it
    counts = pd.DataFrame(parsed).fillna(0).astype(int).sort_index()
    return counts

def main():
    '''
    Command line wrapper for the above
    '''
    args = parse_args()
    seqdf = pd.read_csv(args.seq, index_col=0)
    #stick the tail on if applicable
    if args.tail is not None:
        seqdf.index = [i+args.tail for i in seqdf.index]
    counts = parse_barcodes(args.bam, seqdf)
    #save the output
    counts.to_csv(args.csv)

if __name__ == "__main__":
    main()
