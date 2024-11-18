#!/usr/bin/env python
#PBS -l walltime=47:59:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=64g,vmem=64g
#PBS -m e


"""
Created on Tue Sep 28 12:50:20 2021

@author: Siavash Ghaffari
"""

import pandas as Pd
import numpy as np
import sys
from statistics import mean
from itertools import zip_longest

# Read the command-line argument passed to the interpreter when invoking the script
GENE = sys.argv[1]
Var = sys.argv[2]


# import Stat_applied data
df = Pd.read_csv("P-values_"+GENE+"_NO_Expression_"+Var+"s.csv", sep=',')

# Correcting some columns type
df["Affection_F"] = df["Affection_F"].fillna("[]").apply(lambda x: eval(x))
df["Sex_F"] = df["Sex_F"].fillna("[]").apply(lambda x: eval(x))

# Claening some title destroyed by R usage
df = df.rename(columns={"X..probands":"% probands","X..Unaffected.Sibling":"% Unaffected Sibling","X..1000G":"% 1000G","DGV.Freq":"DGV Freq","DGV.N_Studies":"DGV N_Studies"})

# Fill NA with zeros
if Var=="CNV":
    df["DGV Freq"] = df["DGV Freq"].fillna(0)
    df["DGV N_Studies"] = df["DGV N_Studies"].fillna(0)
elif Var=="SV":
    df["gnomAD_rare_Freq"] = df["gnomAD_rare_Freq"].fillna(0)
    df["gnomAD_common_Freq"] = df["gnomAD_common_Freq"].fillna(0)



# Remove DGV_frequncy less than 0.25
if Var=="CNV":
    df = df[df["DGV Freq"]<0.25]
elif Var=="SV":
    df = df[df["gnomAD_rare_Freq"]<0.25]
    df = df[df["gnomAD_common_Freq"]<1]



# Applying cut-off for all chromosmes except X chromosome
# 1000G test control
df.loc[df["CHROM"]!="chrX","1000_test"] = np.where(df[df["CHROM"]!="chrX"]["N_probands"]>20,df[df["CHROM"]!="chrX"]["pVal_1000G"], 0)
df.loc[df["CHROM"]!="chrX","1000_test"] = np.where(df[df["CHROM"]!="chrX"]["N_1000G"]<=4,0, df[df["CHROM"]!="chrX"]["pVal_1000G"])

# Siblings test Control
df.loc[df["CHROM"]!="chrX","Sibling_test"] = np.where(df[df["CHROM"]!="chrX"]["N_probands"]>20,df[df["CHROM"]!="chrX"]["pVal_siblings"], 0)
df.loc[df["CHROM"]!="chrX","Sibling_test"] = np.where(df[df["CHROM"]!="chrX"]["N_unaffectedsibling"]<=5,0, df[df["CHROM"]!="chrX"]["pVal_siblings"])


# Applying cut-off for X-chromosome
# 1000G test control
df.loc[df["CHROM"]=="chrX","1000_test"] = np.where(df[df["CHROM"]=="chrX"]["N_Control_male"]<=1,0, df[df["CHROM"]=="chrX"]["pVal_1000G"])
df.loc[df["CHROM"]=="chrX","1000_test"] = np.where(df[df["CHROM"]=="chrX"]["N_probands_male"]<=2,df[df["CHROM"]=="chrX"]["pVal_1000G"], 0)

# Siblings test Control
df.loc[df["CHROM"]=="chrX","Sibling_test"] = np.where(df[df["CHROM"]=="chrX"]["N_unaffected_male"]<=2,0, df[df["CHROM"]=="chrX"]["pVal_siblings"])
df.loc[df["CHROM"]=="chrX","Sibling_test"] = np.where(df[df["CHROM"]=="chrX"]["N_probands_male"]<=2,df[df["CHROM"]=="chrX"]["pVal_siblings"], 0)

# Removing very bad data
df = df[df["1000_test"]<0.05].reset_index(drop=True)
df = df[df["Sibling_test"]<0.1].reset_index(drop=True)

# Calculate total P_value  
df["p_value_total"] = np.where(df["CHROM"]=="chrX", np.sqrt(df.pVal_siblings_male*df.pVal_1000G_male),np.sqrt(df.pVal_siblings*df.pVal_1000G))

# Risk Factor evaluation
# for X-Chromosome multiply the effect of affecetion and sex
df["Affection_Xchrom"]=[[
        [1 if AA==BB else 0 for AA, BB in zip_longest(x,y, fillvalue=0)] for x, y in zip_longest(d, c, fillvalue=0)]
        for d, c in zip(df.Affection_F, df.Sex_F)
        ] 

# for all chromosome Find the affection risk score
df["Affection_S"] = df["Affection_F"].apply(lambda x:[mean(a) for a in x])
# X-Chromosome
df["Affection_S_X"] = df["Affection_Xchrom"].apply(lambda x:[mean(a) for a in x])

# Evaluate the final Risk Affection score for all 
df["Affection_Score"] = np.where(df["CHROM"]=="chrX",df["Affection_S_X"].apply(lambda x:np.mean(x)),df["Affection_S"].apply(lambda x:np.mean(x)))

# Write out the final data
df.to_csv("Final_Cleaned_"+GENE+"_"+Var+"s.tsv",sep="\t", index=False)
