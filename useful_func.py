#!/usr/bin/env python
#PBS -l walltime=47:59:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=64g,vmem=64g
#PBS -m e

"""
Created on Tue Sep 28 12:43:33 2021

@author: Siavash Ghaffari
"""

import pandas as Pd

"""
This file consists of 3 useful functions needed to wrap up the final list
"""

# Overlap Algorithm
# Setup the Exons_overlapped function to find each exons overlaps
def Exons_overlapped (Overlap,Exons):
    # Setup a chrom list to simplify the following analysis
    chrom = []
    for i in range(1,23):
        chrom.append ("chr"+str(i))
    chrom.append("chrX")
    chrom.append("chrY")
    lncRNA_list = list(Exons["lncRNA"])    
    A, B, C, D, E, F, G, H, K, L, M, N, O = [],[],[],[],[],[],[],[],[],[],[],[],[]
    for i in range(len(lncRNA_list)):
        lncRNA_gene = Overlap[Overlap["lncRNA_gene"]==lncRNA_list[i]].reset_index(drop=True)
        for j in range(len(lncRNA_gene)):
    
            for k in range(len(Exons["Exons"][i])):
                inter = len(range(max(lncRNA_gene["Start"][j], Exons["Exons"][i][k][0]), min(lncRNA_gene["End"][j], Exons["Exons"][i][k][1])))
                if inter>0:
                #if ((lncRNA_gene["Start"][j] <= Exons["Exons"][i][k][0]) and (lncRNA_gene["End"][j] >= Exons["Exons"][i][k][1])):
                    A.append ((Exons["lncRNA"][i],k+1))
                    B.append (Exons["Exons"][i][k])
                    C.append (inter)
                    D.append(lncRNA_gene["Sample"][j])
                    E.append(lncRNA_gene["SVTYPE"][j])
                    F.append(lncRNA_gene["Start"][j])
                    G.append(lncRNA_gene["End"][j])
                    H.append(lncRNA_gene["Dataset"][j])
                    K.append(lncRNA_gene["CHROM"][j])
                    L.append(lncRNA_gene["High-confidence"][j])
                    M.append(lncRNA_gene["Inheritance"][j])
                    N.append(lncRNA_gene["Description"][j])
                    O.append(lncRNA_gene["SEX"][j])
                                       
                #print (Overlap_chr["Sample"][i],Overlap_chr["lncRNA"][i],Overlap_chr["Description"][i])
                                    
    Exons_overlap = Pd.DataFrame()
    Exons_overlap ["lncRNA,Exon"] = A
    Exons_overlap ["Exon_Co"] = B
    Exons_overlap ["overlap"] = C
    Exons_overlap ["Sample"] = D
    Exons_overlap ["SVTYPE"] = E
    Exons_overlap ["Start"] = F
    Exons_overlap ["End"] = G
    Exons_overlap ["Dataset"] = H
    Exons_overlap ["CHROM"] = K
    Exons_overlap ["High-confidence"] = L
    Exons_overlap ["Inheritance"] = M
    Exons_overlap ["Description"] = N
    Exons_overlap ["SEX"] = O
    
    return Exons_overlap  

# Define a function to calculate the number of samples in each dataset and different CNV inheritance and number high-confidence overlapped CNVs
def Count_adder (df,All):
    dataset = ["MSSNG","SSC"]
    inherit = [ 'P_denovo','Maternal', 'Paternal']
    #All = list(Exons_proband["Inheritance"].unique())
    unknown = [x for x in All if x not in inherit]
    #list(Exons_proband["Inheritance"].unique())
    sex = ["female","male"]
    confid  = ["YES","NO"]
    for i in dataset:
        df["N_"+df.name+"_"+i] = df["Dataset"].apply(lambda x:x.count(i))
    for j in inherit:
        df["N_"+df.name+"_"+j] = df["Inheritance"].apply(lambda x:x.count(j))
    Sum = 0
    for l in unknown:
        Sum = Sum + df["Inheritance"].apply(lambda x:x.count(l))
    df["N_"+df.name+"_Unknown"] = Sum
    for m in sex:
        df["N_"+df.name+"_"+m] = df["SEX"].apply(lambda x:x.count(m))
    for k in confid:
        df["N_"+df.name+"_"+k] = df["High_confidence"].apply(lambda x:x.count(k))
    return df


# Define a function to add diffferent features to our final list
def ADD_DATA (database, main, df, STR):
    database = database.rename(columns={"lncRNA":"lncipediaGeneID"})
    df_database = database.merge(df,
            on=['lncipediaGeneID','Exon_no',"SVTYPE"],
              how='inner')
    df2_database = df_database.drop(columns=["CHROM_x"])
    df3_database = df2_database[["lncipediaGeneID","SVTYPE","Exon_no"]]
    df4_database = df3_database.groupby(["lncipediaGeneID","SVTYPE"]).Exon_no.agg(set)
    df4_database = df4_database.reset_index(drop=False)
    df5_database = df4_database.rename(columns={"Exon_no":STR})
    FINAL = (
    main.merge(df5_database, 
              on=['lncipediaGeneID',"SVTYPE"],
              how='left'))
    return FINAL

# Define a function to calculate the number of male/female and samples in each dataset (just for 1000G and MGRB)
def Count_adder_C (df):
    dataset = ["1000G","MGRB"]

    sex = ["female","male"]
    
    for i in dataset:
        df["N_"+df.name+"_"+i] = df["Dataset"].apply(lambda x:x.count(i))
    
    for m in sex:
        df["N_"+df.name+"_"+m] = df["SEX"].apply(lambda x:x.count(m))
    
    return df

# Define a function to add different features to our final list (for enhnacer_gene intercation)
def ADD_DATA_intercation (database, main, df, gene_lists, STR):
    database = database.rename(columns={"lncRNA":"lncipediaGeneID","Gene_Name":"ASD_Gene"})
    df_database = database.merge(df,
            on=['lncipediaGeneID','Exon_no',"SVTYPE"],
              how='inner')
    df2_database = df_database[df_database["ASD_Gene"].isin(gene_lists)].reset_index(drop=True)
    df2_database["EXON_GENE"] = tuple(zip(df2_database.Exon_no,df2_database.ASD_Gene))
    df3_database = df2_database[["lncipediaGeneID","SVTYPE","EXON_GENE"]]
    df4_database = df3_database.groupby(["lncipediaGeneID","SVTYPE"]).EXON_GENE.agg(set)
    df4_database = df4_database.reset_index(drop=False)
    df5_database = df4_database.rename(columns={"EXON_GENE":STR})
    FINAL = (
    main.merge(df5_database, 
              on=['lncipediaGeneID',"SVTYPE"],
              how='left'))
    return FINAL