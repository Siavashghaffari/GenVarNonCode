#!/usr/bin/env python
#PBS -l walltime=23:59:00
#PBS -l nodes=1:ppn=32
#PBS -l mem=64g,vmem=64g
#PBS -m abe

print ("0")

import pandas as Pd
import numpy as np
from itertools import zip_longest
#import dask.dataframe as dd

"""
THIS CODE finds the overlap between CNVs and lncRNAs from lncipedia for all inheritance
PICK the dataset and inheritance and RUN!
"""
print ("1")
dataset="MSSNG"

CNV_MSSNG = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/SVs.MSSNG.freq_1percent.HQR.tsv", sep='\t', low_memory=False)
CNV_SSC = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/SVs.SSC.freq_1percent.HQR.tsv", sep='\t', low_memory=False)
#print (CNV.columns)
print ("2")

MSSNG_2 = CNV_MSSNG[["sample","chrm","start","end","sv_type","length","method","gene_symbol","gene_egID","gene_symbol_CNVstart","gene_symbol_CNVend","exon_symbol","exon_egID",
        "cds_symbol","cds_egID","gnomAD_pLI","Inheritance","gnomAD_commonSV","gnomAD_rareSV","ASD_gene_lists_gene","ASD_gene_lists_exon","Sex","Relation","Affection","Mother ID",
        "Father ID","Family ID","Predicted ancestry","Family size"]]
SSC_2 = CNV_SSC[["sample","chrm","start","end","sv_type","length","method","gene_symbol","gene_egID","gene_symbol_CNVstart","gene_symbol_CNVend","exon_symbol","exon_egID",
        "cds_symbol","cds_egID","gnomAD_pLI","Inheritance","gnomAD_commonSV","gnomAD_rareSV","ASD_gene_lists_gene","ASD_gene_lists_exon","Sex","Relation","Affection","Mother ID",
        "Father ID","Family ID","Predicted ancestry","Family size"]]



MSSNG_2.to_csv("/hpf/largeprojects/tcagstor/users/sghaffari/SVs.MSSNG.freq_1percent.HQR_Small.tsv", sep='\t')                    
SSC_2.to_csv("/hpf/largeprojects/tcagstor/users/sghaffari/SVs.SSC.freq_1percent.HQR_Small.tsv", sep='\t')

