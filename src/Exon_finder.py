#!/usr/bin/env python
#PBS -l walltime=47:59:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=64g,vmem=64g
#PBS -m e


"""
Created on Fri Aug 20 12:30:51 2021

@author: Siavash Ghaffari
"""

import pandas as Pd
from itertools import zip_longest
import re

def Overlapped_Exons (Overlap):
    """
    This function receives final overlap of lncRNa with CNVs and returns 
    all exons of overlapped lncRNAs genes ((not transcript) with their cordinates
    """
    
    INcRNA = Pd.read_csv("/home/sghaffari/lncipedia_5_2_hg38.bed", sep='\t', header=None)
    
    # Cleaning lncRNA information
    INcRNA.columns= ["CHR", "start", "End", "Transcript_name", "Score", "Strand", "Start_coding", "end_coding", "Color", "blockCount", "blockSizes", "blockStarts"]
    INcRNA["Gene"] = INcRNA["Transcript_name"].str.split(":").map(lambda x: x[0])
    INcRNA["blockSizes"]=INcRNA["blockSizes"].str.rstrip(",")
    INcRNA["blockStarts"]=INcRNA["blockStarts"].str.rstrip(",")
    INcRNA["blocksize_list"] = [list(map(int, i.split(","))) for i in INcRNA["blockSizes"]]
    INcRNA["blockstart_list"] = [list(map(int, i.split(","))) for i in INcRNA["blockStarts"]]
    INcRNA = INcRNA.assign(blockend_list=[
    [x + y for x, y in zip_longest(d, c, fillvalue=0)]
    for d, c in zip(INcRNA.blockstart_list, INcRNA.blocksize_list)
    ])
    
    # Find list of lncRNA that overlapping CNVs
    lncRNA_list = list(Overlap["lncRNA_gene"].unique())
    
    # filter lncRNAs for those ones in the overlaps list
    Final = INcRNA[INcRNA.Gene.isin(lncRNA_list)].reset_index(drop=True)
    
    # Find exons of any lncRNA gene (not transcript) by combining all transcript exons
    H, I = [],[]
    for i in range(len(lncRNA_list)):
        range_i = []
        lncRNA_gene = Final[Final["Gene"]==lncRNA_list[i]].reset_index(drop=True)
        for j in range(len(lncRNA_gene)):
            for k in range(len(lncRNA_gene["blockstart_list"][j])):    
                start = (lncRNA_gene["blockstart_list"][j][k]+lncRNA_gene["start"][j]).astype("int64")
                end = (lncRNA_gene["blockend_list"][j][k]+lncRNA_gene["start"][j]).astype("int64")
                range_i.append((start,end))
        b = []
        for begin,end in sorted(range_i):
            if b and b[-1][1] >= begin - 1:
                b[-1] = (b[-1][0], end)
            else:
                b.append((begin, end))
        H.append(lncRNA_list[i])
        I.append(b)
        
    lncRNA_Exon = Pd.DataFrame()
    lncRNA_Exon ["lncRNA"] = H
    lncRNA_Exon ["Exons"] = I
    return lncRNA_Exon

    