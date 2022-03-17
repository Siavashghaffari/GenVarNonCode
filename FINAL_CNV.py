#!/usr/bin/env python
#PBS -l walltime=47:59:00
#PBS -l nodes=1:ppn=32
#PBS -l mem=64g,vmem=64g
#PBS -m e


"""
Created on Tue Aug 17 12:08:29 2021

@author: Siavash Ghaffari
"""

import pandas as Pd
import numpy as np
from itertools import zip_longest

from CDS_PLI_examiner import CNV_CDS_remover

def Overlap_management(df, Clin=True, CDS=True, Expression=True, CNV_Size_cutoff=3000000):
    """
    Overlap_management function receives intial overlapped CNV and remove non relevant ones step by step
    1- Remove clinically significant CNVs
    2- Remove CNVs also overlapping PCGs (CNV_CDS_remover function)
    3- Remove CNVs overlapping non-Brain Expressed lncRNAs 
    Clin, CDS and Expression are boolean variables specify if we want step 1,2 and 3 in our analysis respectively 
    output is a list of CNVs passed all the filters
    """
    
    # Exclude Inversions with breakpoint outside the gene
    df_2 = df[~((df["Description"]=="Complete overlap") & (df["SVTYPE"]=="INV"))]
    df = df_2.copy()

    # Reading high-quality lncRNAs
    HC = Pd.read_csv("/home/sghaffari/lncipedia_5_2_hc_hg38.bed", sep='\t', header=None)
    HC.columns= ["CHR", "start", "End", "Gene_name", "Score", "Strand", "Start_coding", "end_coding", "Color", "blockCount", "blockSizes", "blockStarts"]
    
    # Create lncRNA gene columns from transcripts (without ":") and size feature from start and end
    df["lncRNA_gene"] = df["lncRNA"].str.split(":").map(lambda x:x[0])
    df["Size"] = df["End"]-df["Start"]+1
    df["Start"] = df["Start"].astype(np.int64)
    df["End"] = df["End"].astype(np.int64)
    
    # Start with Raw Overlaps
    df_new = df.copy()
    
    #Exclude large CNVs
    df_new2 = df_new[df_new.Size<=CNV_Size_cutoff]
    
    # Exclude Clinical significant CNVs
    if Clin==True:
        df2 = Clin_check(df_new2)
    else:
        df2 = df_new2.copy()
    
    # Remove CNVs also overlapping PCGs
    if CDS==True:
        df3 = CNV_CDS_remover (df2, "ASD", PLI_threshold=0.9 , action=True)
    else:
        df3 = df2.copy()
        
    # Group data for lncRNAs instead of smaples to exclude zero exprssed lncRANs
    group = df3.groupby(['lncRNA_gene','SVTYPE'])
    df20 = group.apply(lambda x: Pd.Series({
    "Samples"   :   set(x['Sample']),
    "Transcripts"   :   set(x['lncRNA']),
    "Size"   :   set(x['Size']),
    "CHROM"   :   set(x['CHROM']),
    })
    )
    df20 = df20.reset_index(drop=False)
    
    if Expression==True:
        # Reading Expression Data
        Expression = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/lncRNA_Overlaps_IO/lncRNA_brain_expression_final.tsv", sep='\t')
        Expression_list = Expression["Lincpedia"].to_list()
        df22 = df20[df20["lncRNA_gene"].isin(Expression_list)].reset_index(drop=True)
    else:
        df22 = df20.copy()
        
    df22 = df22.drop_duplicates(subset=['lncRNA_gene','SVTYPE']).reset_index(drop=True)
    
    # Count number of samples for each lncRNA
    df22["N_probands"] = df22["Samples"].str.len()
    
    # Keep only lncRNAs that overlaps iwth at least two samples
    df23= df22[df22["N_probands"]>=1].reset_index(drop=True)
    df24 = df23.sort_values(['N_probands'], ascending=False).reset_index(drop=True)
    df24=df24[["lncRNA_gene","SVTYPE","N_probands"]]
    lncRNA_list = df23["lncRNA_gene"].to_list()
    df4 = df3[df3["lncRNA_gene"].isin(lncRNA_list)].reset_index(drop=True)
    df4 = df4[['Sample', 'lncRNA', 'Description', 'CHROM', 'Start', 'End', 'SVTYPE',
        'lncRNA_gene', 'Size','Inheritance', 'SEX', 'cds_symbol', 'gnomAD_pLI', 'Dataset']]
    
    # label high-confidence lncRNAs
    HC_list = HC["Gene_name"].unique()
    df4["High-confidence"] = np.where(df4["lncRNA"].isin(HC_list),"YES","NO")
    return df4, df24
    

    
def Clin_check(CNV_in):
        #Reading clinical Significant CNVs and cleaning the data
        Genomic_disorder = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/lncRNA_Overlaps_IO/Clinically significant CNVs-MSSNGDB6-SSC-1KG-hg38-20201119_Shared.csv", sep=',',low_memory=False)
        Genomic_disorder_1 = Genomic_disorder[:484]
        Genomic_disorder_1 = Genomic_disorder_1.rename(columns={"sample":"Sample","chr":"CHROM","start":"Start","type":"SVTYPE"})
        Genomic_disorder_1["Start"]=Genomic_disorder_1["Start"].astype(np.int64)
        Genomic_disorder_1["end"]=Genomic_disorder_1["end"].astype(np.int64)
        Genomic_disorder_2 = Genomic_disorder_1[Genomic_disorder_1["Class"].isin(["3-10Mb CNV",'>10Mb CNV','Aneuploidy','Genomic Disorder','Unbalanced Translocation'])].reset_index(drop=True)
        Genomic_disorder = Genomic_disorder_2
        CNV_out = (
        CNV_in.merge(Genomic_disorder, 
        on=['Sample','CHROM','SVTYPE'],
        how='left', 
        indicator=True)
        .query('_merge == "left_only"')
        .drop(columns='_merge')
        )
        CNV_out = CNV_out.reset_index(drop=True)
        CNV_out = CNV_out.rename(columns={"Start_x":"Start"})
        return CNV_out    