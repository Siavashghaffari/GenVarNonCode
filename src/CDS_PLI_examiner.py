#!/usr/bin/env python
#PBS -l walltime=47:59:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=64g,vmem=64g
#PBS -m e

"""
Created on Tue Aug 17 09:30:01 2021

@author: Siavash Ghaffari
"""

import pandas as Pd
from itertools import zip_longest
import re
import numpy as np
from gene_list_pick import Gene_list_pick


def CNV_CDS_remover (CNV, gene_lists, PLI_threshold=0.9 , action=True):
    
    """
    CNV_CDS_remover function removes CNVs also overlapping the PCGs and are believed to impact PCGs more than lncRNA
    (if they are in the ASD (NDD) listand PLI is high enough).
    It removes CNVs in different steps and remaining CNVs are the one and we will continue our analysis with them.
    The input is all CNVs and the output is CNVs that are important for lncRNA analysis
    gene_lists varibale is eithre "ASD" or "NDD" ; PLI_threshold is a parameter to keep CNVs with PLI less than that (the default is 0.9)
    action is a boolean variable tells if should remove theses steps or not (the default is True)
    """
    
    #if the action=True perform the following steps
    if action == True:
        # Replace "-" with NaN to simplfy the analysis later
        CNV["cds_symbol"] = CNV["cds_symbol"].replace("-",np.nan)
        CNV["gnomAD_pLI"] = CNV["gnomAD_pLI"].replace("-",np.nan)       
        


        # check which gene_list we want to remove:ASD or NDD or CP from Gene_list_pick function
        ASD_list = Gene_list_pick (gene_lists)[0]
        
        # Import Barin zero expressed PCGs and make a list out of it
        Brain_zero_expressed_PCG = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/lncRNA_Overlaps_IO/PCGs_Not_Expressed_Brain.tsv", sep="\t",low_memory=False)
        zero_PCG_list = Brain_zero_expressed_PCG["gene"].to_list()
    
        # First Step: Remove CNVS that don't have CDS 
        NO_CDS = CNV[CNV.cds_symbol.isna()]
    
        # Second Step: Remove CNVS that don't have PLI score 
        NO_PLI = CNV[CNV["gnomAD_pLI"].isna()].reset_index(drop=True)
    
        # Third Step: Keep the CNVs that have PLI score
        CNV = CNV.dropna(subset=["gnomAD_pLI","cds_symbol"]).reset_index(drop=True)
        
        # Prepare data for next steps
        #CNV["gnomAD_pLI"] = CNV["gnomAD_pLI"].str.strip("...")
        CNV["PLI_pre"] = CNV["gnomAD_pLI"].str.split("\\|")
        CNV["CDS"] = CNV["cds_symbol"].str.split("\\|")
        
        # Cleaning PLI lists
        temp = "..."
        CNV["EVA"] = CNV["PLI_pre"].apply(lambda x:[temp not in i for i in x])
        CNV = CNV.assign(PLI=[
               [x * y for x, y in zip_longest(d, c, fillvalue=0)]
               for d, c in zip(CNV["PLI_pre"], CNV["EVA"])
                 ])  
        CNV["PLI"] = CNV["PLI"].apply(lambda x:list(filter(None, x)))
        


        # Keep PLIs for only some genes that CNV impcat their CDSs
        E = []
        for index,row in CNV.iterrows(): 
            E.append([i for i in row["PLI"] if i.split(':')[0] in row["CDS"]])
        CNV ["PLI_NOCDS_removed"] = E
    
        # Fourth Step: Remove CNVs that don't have PLI score for their CDS
        NO_CDSPLI = CNV[CNV["PLI_NOCDS_removed"].str.len()==0].reset_index(drop=True)
        # Continue our analysis with therest of CNVs
        CNV = CNV[CNV["PLI_NOCDS_removed"].str.len()>0].reset_index(drop=True)
    
        # Modify our PLI scores regarding to ASD list or being expressed in brain
        A, B, C, D = [],[],[],[]
        for index,row in CNV.iterrows(): 
            A.append([float(i.split(':')[1]) if i.split(':')[1] !="NA" else 0 for i in row["PLI_NOCDS_removed"]])
            B.append([i.split(':')[0] for i in row["PLI_NOCDS_removed"]])
            C.append([1 if i.split(':')[0] in ASD_list else 0 for i in row["PLI_NOCDS_removed"]])
            D.append([0 if i.split(':')[0] in zero_PCG_list else 1 for i in row["PLI_NOCDS_removed"]])
        CNV ["PLI_scores"] = A
        CNV ["PLI_gene"] = B
        CNV ["ASD_list"] = C
        CNV ["zero_PCG_list"] = D          
        CNV["ASD_list_check"] = CNV["ASD_list"].apply(sum)
    
        # Fifth Step: Remove CNVs also overlaps ASD_list from the rest of CNVs
        CNV_new = CNV[CNV["ASD_list_check"]==0]
    
        # Add effect of Barin_Expression on PLI of PCG
        CNV_new["PLI_updated"]=[
        [x * y for x, y in zip_longest(d, c, fillvalue=0)]
        for d, c in zip(CNV_new.PLI_scores, CNV_new.zero_PCG_list)
        ]
        # Find max among all PLI scores
        CNV_new["PLI_scores_max"] = CNV_new.PLI_updated.apply(max)
    
        # Sixth Step:Remove CNVs also overlapping brain expressed PCGs with PLI>0.9
        CNV_new2 = CNV_new[CNV_new["PLI_scores_max"] <= PLI_threshold].reset_index(drop=True)
    
        # Cleaning data
        CNV_new2 = CNV_new2.drop(columns=['PLI', 'PLI_pre','EVA', 'PLI_NOCDS_removed','PLI_scores', 'PLI_gene','ASD_list','zero_PCG_list','CDS','ASD_list_check','PLI_updated','PLI_scores_max'])
        NO_CDSPLI = NO_CDSPLI.drop(columns=['PLI', 'PLI_pre','EVA','CDS','PLI_NOCDS_removed'])
    
        # Union 3 parts of data that we removed before and the rest of CNvs after appropriate filtering
        CNV_out = Pd.concat([NO_CDS, NO_PLI, NO_CDSPLI, CNV_new2]).drop_duplicates().reset_index(drop=True)
    
    else:
        CNV_out = CNV.copy() 
    
    return CNV_out
    