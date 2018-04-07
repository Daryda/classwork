# -*- coding: utf-8 -*-
"""
Created on Sat Apr  7 11:35:10 2018

@author: Darida
"""
from Bio import SeqIO
import argparse

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Needleman–Wunsch algorithm')
    parser.add_argument('-s', '--del_start', help='How many nucleotides cut from start of sequence', metavar="float", type=float, default=0)
    parser.add_argument('-e', '--del_end', help='How many nucleotides cut from end of sequence', metavar="float", type=float, default=0)
    parser.add_argument('-w', '--window', help='How many nucleotides in window', metavar="float", type=float, default=4)
    parser.add_argument('-c', '--cut', help='Delete tresh to..', metavar="float", type=float, default=30)
    parser.add_argument('-i', '--input', help='Input filename for analise', metavar='FILE', required=True)
    parser.add_argument('-p', '--file', help='File for analyse', metavar='FILE', required=True)
    parser.add_argument('-o', '--output', help='Output filename for analise', metavar='FILE', required=True)
    
    args=parser.parse_args()
    
    i = args.input
    start = args.del_start-1
    end = args.del_end-1
    wi=args.window
    cut=args.cut
    o = args.output
    
    
    records = (rec[start:end] for rec in SeqIO.parse(open(i), "fastq"))
    handle = open("p.fastq", "w")
    count = SeqIO.write(records, handle, "fastq")
    handle.close()
## Open FASTA files
#    with open(i, "rU") as handle:
#        seq_list=[]   
#        for record in SeqIO.parse(handle, "fastq"):
#            seq_list.append(str(record.seq))
#    nseq=len(seq_list)
#    seq_st=seq_list.remove(seq_list[:start-1])
#    seq=seq_st.remove(seq_st[end+1:])

import matplotlib.pyplot as plt

    count = 0
    for rec in SeqIO.parse("p.fastq", "fastq"):
        count += 1
    print "%i reads" % count
    good_reads = (rec for rec in \
                  SeqIO.parse("p.fastq", "fastq") \
                  if min(rec.letter_annotations["phred_quality"]) >= cut)
    count = SeqIO.write(good_reads, "good_quality.fastq", "fastq")

    print "Saved %i reads" % count
plt.bar(list(dictionary.keys()), dictionary.values(), color='g')
plt.show()

