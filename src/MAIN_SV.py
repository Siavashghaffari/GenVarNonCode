#!/usr/bin/env python
#PBS -l walltime=47:59:00
#PBS -l nodes=1:ppn=32
#PBS -l mem=64g,vmem=64g
#PBS -m e

"""
Created on Thu Aug 19 10:51:51 2021

@author: Siavash Ghaffari
"""

import pandas as Pd
from itertools import zip_longest
import re
import math
import numpy as np
import scipy.stats as stats

from FINAL_CNV import Overlap_management
from Exon_finder import Overlapped_Exons
from features import Features
from gene_list_pick import Gene_list_pick
from useful_func import *
from add_annotations_SV import ADD_annotations

"""
    Overlap_management function receives intial overlapped CNV and remove non relevant ones step by step
    1- Remove clinically significant CNVs
    2- Remove CNVs also overlapping PCGs (CNV_CDS_remover function)
    3- Remove CNVs overlapping non-Brain Expressed lncRNAs 
    Clin, CDS and Expression are boolean variables specify if we want step 1,2 and 3 in our analysis respectively 
    output is a list of CNVs passed all the filters
"""

# Set parameters here and run
clin = False
cds = True
expression = False
CNV_size_cutoff = 3000000

# What to label print out file
if (clin, cds, expression) == (True, True, True):
    case = ""
elif (clin, cds, expression) == (False, False, False):
    case = "_NO_CLINICAL_CDS_EXPRESSION"
elif (clin, cds, expression) == (False, True, False):
    case = "_NO_EXPRESSION"
elif (clin, cds, expression) == (False, False, True):
    case = "_NO_CDS"
elif (clin, cds, expression) == (False, True, True):
    case = "_NO_CLINICAL"
else:
    case = ""

# Objectify and get all probands and unaffected siblings overlaps, as well as 1000G
MSSNG = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/lncRNA_Overlaps_IO/lncRNA_Overlapped_SV_MSSNG_.tsv", sep='\t',low_memory=False)
SSC = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/lncRNA_Overlaps_IO/lncRNA_Overlapped_SV_SSC_.tsv", sep='\t',low_memory=False)
MSSNG_unaffected = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/lncRNA_Overlaps_IO/lncRNA_Overlapped_SV_MSSNG_unaffected.tsv", sep='\t',low_memory=False)
SSC_unaffected = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/lncRNA_Overlaps_IO/lncRNA_Overlapped_SV_SSC_unaffected.tsv", sep='\t',low_memory=False)
#import 1000G as control
CNV_1kG = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/lncRNA_Overlaps_IO/lncRNA_Overlapped_SV_1000G_.tsv", sep='\t',low_memory=False)
#CNV_MGRB = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/lncRNA_Overlaps_IO/lncRNA_Overlapped_MGRB_.tsv", sep='\t',low_memory=False)
# import affected parents
MSSNG_parents = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/lncRNA_Overlaps_IO/lncRNA_Overlapped_SV_MSSNG_parents.tsv", sep='\t',low_memory=False)

# Specify datasets and combine data come from different database
SSC ["Dataset"] = "SSC"
MSSNG ["Dataset"] = "MSSNG"
df_case = Pd.concat([MSSNG,SSC], ignore_index=True)
SSC_unaffected ["Dataset"] = "SSC"
MSSNG_unaffected ["Dataset"] = "MSSNG"
df_unaffected = Pd.concat([MSSNG_unaffected,SSC_unaffected], ignore_index=True)
CNV_1kG ["Dataset"] = "1000G"
#CNV_MGRB ["Dataset"] = "MGRB"
#CNV_1000G = Pd.concat([CNV_1kG,CNV_MGRB], ignore_index=True)
MSSNG_parents["Dataset"] = "MSSNG"

# Replace "-" with NaN to simplfy the analysis later
#df_case["cds_symbol"] = df_case["cds_symbol"].replace('-', np.nan)
#df_unaffected["cds_symbol"] = df_unaffected["cds_symbol"].replace('-', np.nan)
#CNV_1kG["cds_symbol"] = CNV_1kG["cds_symbol"].replace('-', np.nan)
#MSSNG_parents["cds_symbol"] = MSSNG_parents["cds_symbol"].replace('-', np.nan)

# Calling overlap_management function for all probands, unaffected and 1000G to get the overlaps and counts
Overlap_probands, Proband_counter = Overlap_management(df_case, clin, cds, expression, CNV_size_cutoff)
Overlap_unaffected, unaffected_counter = Overlap_management(df_unaffected, clin, cds, expression, CNV_size_cutoff)
Overlap_1000G, C1000G_counter = Overlap_management(CNV_1kG, clin, cds, expression, CNV_size_cutoff)
Overlap_parents, parents_counter = Overlap_management(MSSNG_parents, clin, cds, expression, CNV_size_cutoff)

# Replace Na with "-" to simplfy the analysis later
Overlap_probands = Overlap_probands.replace(np.nan, '-', regex=True)
Overlap_unaffected = Overlap_unaffected.replace(np.nan, '-', regex=True)
Overlap_parents = Overlap_parents.replace(np.nan, '-', regex=True)
# Calling Exon_finder funstion to get all the exons info
Exons = Overlapped_Exons(Overlap_probands)

# Call the function for probands
Exons_proband  = Exons_overlapped (Overlap_probands,Exons)
# Drop duplicates from the overlap database
Exons_proband = Exons_proband.drop_duplicates().reset_index(drop=True)
# Get a copy of overlap info for future analysis
Overlap_info = Exons_proband.copy()
# Drop the overlp description column to clean the data
Exons_proband = Exons_proband.drop(columns=["Description"])
Exons_proband = Exons_proband.drop_duplicates().reset_index(drop=True)
# Drop "NO_high-confidence" overlap while the high-confidence one is available for the same overlap
Exons_proband_2 = Exons_proband[(Exons_proband.duplicated(['lncRNA,Exon','Sample','SVTYPE',"Dataset"],keep=False)) & (Exons_proband["High-confidence"]=="NO")]
# Remove Exons_proband_2 form Exons_proband
Exons_proband_3 = Pd.concat([Exons_proband, Exons_proband_2, Exons_proband_2]).drop_duplicates(keep=False)
# Save the final and cleaned data in the final variable
Exons_proband = Exons_proband_3.copy()
# Call the function for unaffected siblings
Exons_unaffected  = Exons_overlapped (Overlap_unaffected,Exons)
# Drop duplicated for unaffcted siblings same as what we did for probands
Exons_unaffected = Exons_unaffected.drop_duplicates().reset_index(drop=True)
# Get a copy from overlap information of unaffcted siblings
Overlap_info_unaffected = Exons_unaffected.copy()
# Drop the overlap description coulmn from unaffected siblings
Exons_unaffected = Exons_unaffected.drop(columns=["Description"])
Exons_unaffected = Exons_unaffected.drop_duplicates().reset_index(drop=True)
Exons_unaffected_2 = Exons_unaffected[(Exons_unaffected.duplicated(['lncRNA,Exon','Sample','SVTYPE',"Dataset"],keep=False)) & (Exons_unaffected["High-confidence"]=="NO")]
Exons_unaffected_3 = Pd.concat([Exons_unaffected, Exons_unaffected_2, Exons_unaffected_2]).drop_duplicates(keep=False)
# Save the final and cleaned data in the final variable
Exons_unaffected = Exons_unaffected_3.copy()
# Call the function for 1000G 
Exons_1000G  = Exons_overlapped (Overlap_1000G,Exons)
# Drop duplicated for 1000G same as what we did for probands
Exons_1000G = Exons_1000G.drop_duplicates().reset_index(drop=True)
# Drop the overlap description coulmn from 1000G
Exons_1000G = Exons_1000G.drop(columns=["Description"])
Exons_1000G = Exons_1000G.drop_duplicates().reset_index(drop=True)
Exons_1000G_2 = Exons_1000G[(Exons_1000G.duplicated(['lncRNA,Exon','Sample','SVTYPE'],keep=False)) & (Exons_1000G["High-confidence"]=="NO")]
Exons_1000G_3 = Pd.concat([Exons_1000G, Exons_1000G_2, Exons_1000G_2]).drop_duplicates(keep=False)
Exons_1000G = Exons_1000G_3.copy()

# Cleaning the exon column for all databases
Exons_proband["Exon"] = Exons_proband['lncRNA,Exon'].map(lambda x:str(x))
Exons_unaffected["Exon"] = Exons_unaffected['lncRNA,Exon'].map(lambda x:str(x))
Exons_1000G["Exon"] = Exons_1000G['lncRNA,Exon'].map(lambda x:str(x))

# Probands overlaps grouped by lncRNA, exon and SVTYPE
group = Exons_proband.groupby(["Exon","SVTYPE"])
df = group.apply(lambda x: Pd.Series({
    "Samples"   :   set(x['Sample']),
    "Exon_Co"   :  x['Exon_Co'].unique()[0],
    "CHROM"   :  x['CHROM'].unique()[0],
    "Dataset"   :   list(x['Dataset']),
    "Inheritance" : list(x['Inheritance']),
    "High_confidence" : list(x["High-confidence"]),
    "SEX" : list(x["SEX"])
})
)
df = df.reset_index()
df["N_probands"] = df["Samples"].str.len()
df2= df[df["N_probands"]>=2].reset_index(drop=True)
Probands = df2.sort_values(['N_probands'], ascending=False).reset_index(drop=True)

# Call the Count_adder function for probands
All = list(Exons_proband["Inheritance"].unique())
Probands.name = "probands"
Probands = Count_adder(Probands,All)

# Unaffected siblings overlaps grouped by lncRNA, exon and SVTYPE
group2 = Exons_unaffected.groupby(["Exon","SVTYPE"])
df3 = group2.apply(lambda x: Pd.Series({
    "Samples"   :   set(x['Sample']),
    "Exon_Co"   :  x['Exon_Co'].unique()[0],
    "CHROM"   :  x['CHROM'].unique()[0],
    "Dataset"   :   list(x['Dataset']),
    "Inheritance" : list(x['Inheritance']),
    "High_confidence" : list(x["High-confidence"]),
    "SEX" : list(x["SEX"])
})
)
df3 = df3.reset_index()
df3["N_probands"] = df3["Samples"].str.len()
df4= df3[df3["N_probands"]>=1].reset_index(drop=True)
unaffected = df4.sort_values(['N_probands'], ascending=False).reset_index(drop=True)

# Call the Count_adder function for unaffected siblings
unaffected.name = "unaffected"
unaffected = Count_adder(unaffected,All)

# 1000G overlaps grouped by lncRNA, exon and SVTYPE
group = Exons_1000G.groupby(["Exon","SVTYPE"])
df5 = group.apply(lambda x: Pd.Series({
    "Samples"   :   set(x['Sample']),
    "Exon_Co"   :  x['Exon_Co'].unique()[0],
    "CHROM"   :  x['CHROM'].unique()[0],
    "Dataset"   :   list(x['Dataset']),
    "SEX" : list(x["SEX"]),
    "High_confidence" : list(x["High-confidence"])
})
)
df5 = df5.reset_index()
df5["N_probands"] = df5["Samples"].str.len()
df6 = df5[df5["N_probands"]>=1].reset_index(drop=True)
C_1000G = df6.sort_values(['N_probands'], ascending=False).reset_index(drop=True)

# Call the Count_adder function for 1000G
C_1000G.name = "Control"
C_1000G = Count_adder_C(C_1000G)

# parents overlaps grouped by lncRNA, exon and SVTYPE
group = Overlap_parents.groupby(["lncRNA_gene","SVTYPE"])
df7 = group.apply(lambda x: Pd.Series({
    "Samples_parents"   :   set(x['Sample'])
})
)
df7 = df7.reset_index()
df7 = df7.rename(columns={"lncRNA_gene":"lncRNA"})

# Merge probands and unaffcted siblings databases
Final3 = Probands.merge(unaffected,
            on=['Exon', 'SVTYPE'],
              how='left')
Final2 = Final3.merge(C_1000G,
            on=['Exon', 'SVTYPE'],
              how='left')


# Cleaning data
Final2["Exon"] = Final2["Exon"].str.rstrip(")")
Final2["Exon"] = Final2["Exon"].str.lstrip("(")
Final2[['lncRNA', 'Exon_no']] = Final2['Exon'].str.split(',', expand=True)
Final2["lncRNA"] = Final2["lncRNA"].str.lstrip(" \'")
Final2["lncRNA"] = Final2["lncRNA"].str.rstrip(" \'")

# ADd parents datasets
Final = Final2.merge(df7,
            on=['lncRNA', 'SVTYPE'],
              how='left')

# Cleaning Data
Final["N_probands_y"] = Final["N_probands_y"].fillna(0)
Final["N_probands_x"] = Final["N_probands_x"].fillna(0)
Final["N_probands"] = Final["N_probands"].fillna(0)
Final = Final.drop(columns=["Exon_Co","CHROM","High_confidence"])
Final = Final.rename(columns={"Samples_x":"Samples_probands","Samples_y":"Samples_unaffectedsibling","Exon_Co_x":"Exon_Co", "overlap_x":"overlap","Dataset_x":"Dataset_probands",
                             "Dataset_y":"Dataset_unaffectedsibling","N_probands_x":"N_probands","N_probands_y":"N_unaffectedsibling","CHROM_x":"CHROM",
                             "N_probands_YES":"N_probands_highConfidence","N_unaffected_YES":"N_unaffected_highConfidence",
                             "Samples":"Samples_1000G","N_probands":"N_1000G"})
Final["% probands"] = Final["N_probands"]/6858*100
Final["% Unaffected Sibling"] = Final["N_unaffectedsibling"]/2175*100
Final["% 1000G"] = Final["N_1000G"]/(2504+4010)*100
Final = Final[['lncRNA', 'Exon_no',  'Exon_Co','SVTYPE', 'CHROM', 'N_probands','N_unaffectedsibling','N_1000G','% probands','% Unaffected Sibling','% 1000G',
        'Samples_probands', 'N_probands_MSSNG', 'N_probands_SSC', 'N_probands_P_denovo', 'N_probands_Maternal', 'N_probands_Paternal','N_probands_Unknown',
        'N_probands_female', 'N_probands_male', 'N_probands_highConfidence', 'Samples_unaffectedsibling','N_unaffected_MSSNG', 'N_unaffected_SSC', 'N_unaffected_P_denovo',
        'N_unaffected_Maternal','N_unaffected_Paternal','N_unaffected_Unknown','N_unaffected_female','N_unaffected_male','N_unaffected_highConfidence','Samples_1000G', 'N_Control_female','N_Control_male', 'N_Control_1000G', 'N_Control_MGRB','Samples_parents']]

# Add column of the whole gene overlap not only each exons
Proband_counter = Proband_counter.rename(columns={"lncRNA_gene":"lncRNA","N_probands":"N_probands_genes"})
unaffected_counter = unaffected_counter.rename(columns={"lncRNA_gene":"lncRNA","N_probands":"N_unaffected_genes"})
C1000G_counter = C1000G_counter.rename(columns={"lncRNA_gene":"lncRNA","N_probands":"N_1000G_genes"})
Final_new = (
    Final.merge(Proband_counter, 
              on=['lncRNA','SVTYPE'],
              how="left")
)
Final_new2 = (
    Final_new.merge(unaffected_counter, 
              on=['lncRNA','SVTYPE'],
              how="left")
)
Final = (
    Final_new2.merge(C1000G_counter, 
              on=['lncRNA','SVTYPE'],
              how="left")
)
Final["N_probands_genes"] = Final["N_probands_genes"].fillna(0)
Final["N_unaffected_genes"] = Final["N_unaffected_genes"].fillna(0)
Final["N_1000G_genes"] = Final["N_1000G_genes"].fillna(0)

# Add data of ensemble id or alias name
Final = Final.rename(columns={"lncRNA":"lncipediaGeneID"})
Expression = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/lncRNA_Overlaps_IO/lncRNA_brain_expression_final.tsv", sep='\t')
Expression = Expression.rename(columns={"Lincpedia":"lncipediaGeneID"})
Final_2 = (
    Final.merge(Expression, 
              on=['lncipediaGeneID'],
              how="left")
)
Ensemble = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/lncRNA_Overlaps_IO/lncipedia_5_2_ensembl_92_genes.txt", sep='\t')
Ensemble["ensemblGeneID"] = Ensemble["ensemblGeneID"].str.split(".").map(lambda x:x[0])
GTEX = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/lncRNA_Overlaps_IO/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", sep='\t',skiprows=2)
GTEX["ensemblGeneID"] = GTEX["Name"].str.split(".").map(lambda x:x[0])
GTEX = GTEX[["ensemblGeneID","Description"]]
GTEX = GTEX.rename(columns={"Description":"GENE_symbol"})
mart = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/lncRNA_Overlaps_IO/mart_export.txt", sep='\t')
mart = mart.drop_duplicates().reset_index(drop=True)
mart = mart.rename(columns={"Gene stable ID":"ensemblGeneID"})
Convertor_2 = (
    Ensemble.merge(GTEX, 
              on=['ensemblGeneID'],
              how='left'))
Convertor = (
    Convertor_2.merge(mart, 
              on=['ensemblGeneID'],
              how='left'))
Convertor = Convertor.drop(columns=["RefSeq ncRNA ID"])
Convertor = Convertor.drop_duplicates().reset_index(drop=True)
Convertor_list = Convertor["lncipediaGeneID"].unique()

C,D,E = [],[],[]
for i in range(len(Convertor_list)):
    Convertor_lncRNA = Convertor[Convertor["lncipediaGeneID"]==Convertor_list[i]].reset_index(drop=True)
    A,B = [],[]
    for j in range(len(Convertor_lncRNA)):
        A.append(Convertor_lncRNA["ensemblGeneID"][j])
        B.append(Convertor_lncRNA["GENE_symbol"][j])
        B.append(Convertor_lncRNA["Gene name"][j])
    C.append(Convertor_list[i])
    D.append(A)
    B = [str(x) for x in B]
    B = [x for x in B if x != 'nan']
    B = np.asarray(B)
    E.append(np.unique(B))
Convertor_new = Pd.DataFrame()
Convertor_new ["lncipediaGeneID"] = C
Convertor_new ["ensemblGeneID"] = D
Convertor_new ["Gene_Name"] = E

Final_3 = (
    Final_2.merge(Convertor_new, 
              on=['lncipediaGeneID'],
              how='left'))
Final_3 = Final_3[['lncipediaGeneID', 'ensemblGeneID', 'Gene_Name','Exon_no','Exon_Co','SVTYPE', 'CHROM', 'N_probands','N_unaffectedsibling','N_1000G','% probands',
                   '% Unaffected Sibling','% 1000G','Samples_probands','N_probands_MSSNG', 'N_probands_SSC', 'N_probands_P_denovo','N_probands_Maternal', 'N_probands_Paternal',
                   'N_probands_Unknown','N_probands_female', 'N_probands_male', 'N_probands_highConfidence','Samples_unaffectedsibling','N_unaffected_MSSNG', 'N_unaffected_SSC',
                   'N_unaffected_P_denovo', 'N_unaffected_Maternal','N_unaffected_Paternal', 'N_unaffected_Unknown', 'N_unaffected_female','N_unaffected_male','N_unaffected_highConfidence',
                   'Samples_1000G','N_Control_female','N_Control_male', 'N_Control_1000G', 'N_Control_MGRB','Expression','N_probands_genes', 'N_unaffected_genes', 'N_1000G_genes','Samples_parents']]
# Complete overlap investigation
Overlap_info2 = Overlap_info[Overlap_info["High-confidence"]=="YES"]
Overlap_info2 = Overlap_info2[Overlap_info2["Description"]=="Complete overlap"]
Overlap_info2["Exon"] = Overlap_info2['lncRNA,Exon'].apply(lambda x:str(x))
Overlap_info2["Exon"] = Overlap_info2["Exon"].str.rstrip(")")
Overlap_info2["Exon"] = Overlap_info2["Exon"].str.lstrip("(")
Overlap_info2[['lncipediaGeneID', 'Exon_no']] = Overlap_info2['Exon'].str.split(',', expand=True)
Overlap_info2["lncipediaGeneID"] = Overlap_info2["lncipediaGeneID"].str.lstrip(" \'")
Overlap_info2["lncipediaGeneID"] = Overlap_info2["lncipediaGeneID"].str.rstrip(" \'")
group = Overlap_info2.groupby(["lncipediaGeneID","SVTYPE"])
df20 = group.apply(lambda x: Pd.Series({
    "Samples"   :   set(x['Sample']),
})
)
df20 = df20.reset_index(drop=False)
df20["N_Complete_overlap"] = df20["Samples"].str.len()
Final_4 = (
    Final_3.merge(df20, 
              on=['lncipediaGeneID','SVTYPE'],
              how='left'))
Final_4 = Final_4.sort_values(["N_probands"],ascending=False)
Final_4 = Final_4[['lncipediaGeneID', 'ensemblGeneID', 'Gene_Name','Exon_no',  'Exon_Co','SVTYPE', 'CHROM', 'N_probands','N_unaffectedsibling','N_1000G','% probands','% Unaffected Sibling','% 1000G',
         'Samples_probands','N_probands_MSSNG', 'N_probands_SSC', 'N_probands_P_denovo','N_probands_Maternal', 'N_probands_Paternal', 'N_probands_Unknown','N_probands_female', 'N_probands_male',
        'N_probands_highConfidence','N_Complete_overlap','Samples_unaffectedsibling','N_unaffected_MSSNG', 'N_unaffected_SSC','N_unaffected_P_denovo', 'N_unaffected_Maternal','N_unaffected_Paternal', 
        'N_unaffected_Unknown', 'N_unaffected_female','N_unaffected_male','N_unaffected_highConfidence','Samples_1000G', 'N_Control_female','N_Control_male', 'N_Control_1000G', 'N_Control_MGRB','N_probands_genes', 'N_unaffected_genes',
        'N_1000G_genes','Samples_parents','Expression']]

print ("0")

# Export data (still feature columns not added)
Final_4.to_csv("Exons_probands+unaffected_SV"+case+".tsv",sep="\t",index=False)


# Cleaning the list for publishing
Final_5 = Final_4.drop(columns=["Exon_Co"])
#Final_5 = Final_5.sort_values([ 'pVal_1000G' , 'pVal_siblings'], ascending=[1,1]).reset_index(drop=True)
Exons_final = Final_5.groupby(["lncipediaGeneID","SVTYPE"]).Exon_no.agg(set)
Exons_final = Exons_final.reset_index(drop=False)

print ("1")

Final_6 = (
    Final_5.merge(Exons_final, 
              on=['lncipediaGeneID',"SVTYPE"],
              how='left'))
Final_7 = Final_6.drop_duplicates(["lncipediaGeneID","SVTYPE"]).reset_index(drop=True)
Final_7 = Final_7.rename(columns={"Exon_no_y":"Exon_no"})
Final_7 = Final_7.drop(columns=["Exon_no_x"])

print("before Conservation")
# Add different features to our list to evaluate importance
# Import different conserved elements 
GERP = Pd.read_csv("/home/sghaffari/Conserved_elements_hg38.tsv", sep='\t')
VISTA = Pd.read_csv("/home/sghaffari/VISTA_Enhancer_hg38.tsv", sep='\t')
Ultra = Pd.read_csv("/home/sghaffari/Ultra_conserved_hg38.tsv", sep='\t')
# Create a Features instance

features = Features(Final_4)

print("after feature call")

# Call the Conserved_exons method for features instance to create GERP conserved elemnets feature
conserved_exons = features.Conserved_exons(GERP)
# Add the final conserved elements as a column to our database
Final_8 = ADD_DATA(conserved_exons, Final_7, Final_4, "Conserved_exon")

print ("after ADD DATA")

# call the Conserved_exons method for features instance to create Ultra conserved elemnets feature
ultra_exons = features.Conserved_exons(Ultra)
# Add the final Ultra elements as a column to our database
Final_9 = ADD_DATA(ultra_exons, Final_8, Final_4, "Ultra_exons")

# call the Conserved_exons method for features instance to create VISTA conserved elemnets feature
vista_exons = features.Conserved_exons(VISTA, if_Enhancer=True)
# Add the final VISTA elements as a column to our database
Final_10 = ADD_DATA(vista_exons, Final_9, Final_4, "Vista_exons")

print("Conservation Done")

# Import different Enhancer Atlas elements
Astrocyte = Pd.read_csv("/home/sghaffari/Enhancer_Atlas_Astrocyte_hg38.tsv",sep="\t")
Cerebellum = Pd.read_csv("/home/sghaffari/Enhancer_Atlas_cerebellum_hg38.tsv",sep="\t")
ESC = Pd.read_csv("/home/sghaffari/Enhancer_Atlas_ESC_hg38.tsv",sep="\t")
Fetal = Pd.read_csv("/home/sghaffari/Enhancer_Atlas_fetal_hg38.tsv",sep="\t")
Astrocyte_EP = Pd.read_csv("/home/sghaffari/Astrocyte_EP_Cleaned_hg38.tsv",sep="\t")
ESC_EP = Pd.read_csv("/home/sghaffari/ESC_neuron_EP_Cleaned_hg38.tsv",sep="\t")
Cerebellum_EP = Pd.read_csv("/home/sghaffari/Cerebellum_EP_Cleaned_hg38.tsv",sep="\t")

# Create an automated database adder for different enhancers
ENHANCERS = [Astrocyte, Cerebellum, ESC, Fetal]
Astrocyte.name  = "Astrocyte"
Cerebellum.name  = "Cerebellum"
ESC.name  = "ESC"
Fetal.name  = "Fetal"

# Add the final different enhancer elements as a column to our database
df_enh = Final_10.copy()
for enh in ENHANCERS:
    enhancer_exons = features.Conserved_exons(enh, if_Score=True)
    df_enh = ADD_DATA(enhancer_exons, df_enh, Final_4, enh.name+"_exons")
Final_11 = df_enh.copy()
    
# Create an automated database adder for different enhancers_gene interaction
ENHANCERS_intercation = [Astrocyte_EP, ESC_EP, Cerebellum_EP]
Astrocyte_EP.name  = "Astrocyte_EP"
ESC_EP.name  = "ESC_EP"
Cerebellum_EP.name  = "Cerebellum_EP"

# check in both ASD or NDD Gene_lists
GENE_LISTS = ["ASD","NDD"]

# Add the final different enhancer elements as a column to our database
df_enh = Final_11.copy()
for enh in ENHANCERS_intercation:
    for gl in GENE_LISTS:
        enhancer_exons = features.Conserved_exons(enh, if_intercation=True)
        # check which gene_list we want to remove:ASD or NDD or CP from Gene_list_pick function
        ASD_list = Gene_list_pick (gl)[0]
        df_enh = ADD_DATA_intercation(enhancer_exons, df_enh, Final_4, ASD_list, enh.name+"_"+gl)
Final_12 = df_enh.copy()


print("Enhancer Done")

# CO_EXPRESSION NETWORK USING LncBook and lncRNAKB DATABASE
# Import lncRNA Conversion for LncBook and lncRNAKB
LncBook = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/lncRNA_Overlaps_IO/LNCpedia_LncBook_conversion.tsv",sep="\t")
lncRNAKB = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/lncRNA_Overlaps_IO/LNCpedia_lncRNAKB_conversion.tsv",sep="\t")
GTEX = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/lncRNA_Overlaps_IO/GTEX_convertor.tsv", sep='\t')

# Call the Make_Conversion_list_for_Coexp method for features instance to create LncBook and lncRNAKB converted database 
LncBook_Final_list = features.Make_Conversion_list_for_Coexp(LncBook,lncRNAKB,GTEX)

# LncBook CO_EXPRESSION
CELL_LINE = ["cerebellum", "brain","HPA_pcc"]
df_1 = Pd.DataFrame() 
for cl in CELL_LINE:
    coexp = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/LncBook_Coexp/"+cl+"_select_out.txt", sep='\t',low_memory=False)
    coexp_exons_LncBook = features.Co_exp_LncBook(coexp, "ASD")
    coexp_exons_LncBook["Dataset"] = cl
    df_1 = Pd.concat([df_1,coexp_exons_LncBook]).reset_index(drop=True)
df_2 = df_1[["Lincpedia","gene","r","p_value"]]
df_2["r"] = df_2["r"].apply(lambda x:abs(x))
df_2 = df_2.rename(columns={"Lincpedia":"lncipediaGeneID","gene":"ASD_gene"})
group = df_2.groupby('lncipediaGeneID')
df_3 = group.apply(lambda x: Pd.Series({
    "ASD Gene"   :   list(x['ASD_gene']),
    "r"   :   list(x['r']),
    "p_value" : list(x['p_value'])
})
)
df_4=df_3.reset_index(drop=False)
df_4["Connectivity"] = df_4["r"].apply(lambda x:sum(x))
df_4["N_ASDgenes"] = df_4["ASD Gene"].str.len()
df_4["Normalized_Connectivity"] = (df_4["Connectivity"]-min(df_4["Connectivity"]))/(max(df_4["Connectivity"])-min(df_4["Connectivity"]))
df_4 = df_4[["lncipediaGeneID","N_ASDgenes","Normalized_Connectivity"]]

# Merge LncBook Co-Exp data with last Final data to make new columns for co-expression data
Final_13 = (
    Final_12.merge(df_4, 
              on=['lncipediaGeneID'],
              how='left'))

# lncrnakb CO_EXPRESSION
coexp_lncrnakb = features.Co_exp_lncrnakb("ASD")
df_lncrnakb = coexp_lncrnakb[["Lincpedia","ASD_gene","corr"]]
df_lncrnakb_1 = df_lncrnakb.rename(columns={"Lincpedia":"lncipediaGeneID"})

group = df_lncrnakb_1.groupby('lncipediaGeneID')
df_lncrnakb_2 = group.apply(lambda x: Pd.Series({
    "ASD Gene"   :   set(x['ASD_gene']),
    "corr"   :   list(x['corr']),
})
)
df_lncrnakb_3 = df_lncrnakb_2.reset_index(drop=False)
df_lncrnakb_3["Connectivity"] = df_lncrnakb_3["corr"].apply(lambda x:sum(x))
df_lncrnakb_3["N_ASDgenes"] = df_lncrnakb_3["ASD Gene"].str.len()
df_lncrnakb_3["Normalized_Connectivity"] = (df_lncrnakb_3["Connectivity"]-min(df_lncrnakb_3["Connectivity"]))/(max(df_lncrnakb_3["Connectivity"])-min(df_lncrnakb_3["Connectivity"]))
df_lncrnakb_4 = df_lncrnakb_3[["lncipediaGeneID","ASD Gene","corr","N_ASDgenes","Normalized_Connectivity"]]

# Merge lncrnakb Co-Exp data with last Final data to make new columns for the adult co-expression data
Final_14 = (
    Final_13.merge(df_lncrnakb_4, 
              on=['lncipediaGeneID'],
              how='left'))

print ("before add_annotation")


# Create a file with lncRNA annotations in the CNV file
Final_15 = ADD_annotations (Final_14, Overlap_probands, Overlap_unaffected, Overlap_parents, case)

# Calculate p_values for three groups of controls
#Final_15["N_DGV"] = Final_15["DGV Freq"]*100
#Final_15["all_DGV"] = Final_15["N_DGV"].sum()

Final_15["all_probands"] = 6858
Final_15["all_siblings"] = 2175
Final_15["all_1000G"] = 2504

Final_15["male_probands"] = 2094+3533
Final_15["male_siblings"] = 928+105
Final_15["male_1000G"] = 1233

Final_15.to_csv("lncRNA_Overlapping_SVs"+case+".tsv",sep="\t",index=False)


# p_value for siblings control
#Final_15['pVal_siblings'] = Final_15.apply(
#    lambda r: (stats.boschloo_exact([[r.all_probands,r.N_probands], [r.all_siblings, r.N_unaffectedsibling]], alternative="less")).pvalue,
#    axis=1)
# p_value for 1000G control
#Final_15['pVal_1000G'] = Final_15.apply(
#    lambda r: (stats.boschloo_exact([[r.all_probands,r.N_probands], [r.all_1000G, r.N_1000G]], alternative="less")).pvalue,
#    axis=1)
# p_value for DGV control
#Final_15['pVal_DGV'] = Final_15.apply(
#    lambda r: (stats.boschloo_exact([[r.all_probands,r.N_probands], [r.all_DGV, r.N_DGV]], alternative="less")).pvalue,
#    axis=1)

#Final_15.to_csv("lncRNA_Overlapping_SVs_with_pvalue"+case+".tsv",sep="\t",index=False)









