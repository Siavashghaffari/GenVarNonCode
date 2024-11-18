#!/usr/bin/env python
#PBS -l walltime=47:59:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=64g,vmem=64g
#PBS -m e

"""
Created on Thu Sep  9 11:17:23 2021

@author: Siavash Ghaffari
"""

import pandas as Pd
import numpy as np

def Gene_list_pick (gene_lists):
    """
    Gene_list_pick function Assigns our gene_list to one of "ASD", "NDD", "CP" 
    gene_lists is a string variable which can be one of above-mentioned variables, it returns the list according to the input avariable
    None as an input will return an ampty list
    if the variable is not defined as one of those three and None it will give an error
    it returns: 
    ASD_list : gene name list
    Ensemble_list : ensemble ID of our gene list
    DF : a dataframe that converts gene name to ensemble ID 
    """
    if gene_lists == "ASD":
        # reading Final ASD (TADA) genes and make a list out of it
        ASD = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/Gene_lists/ASD_updated.tsv", sep="\t",low_memory=False)
        ASD_list = list(ASD["gene"].unique())
        Ensemble_list = list(ASD["Gene stable ID"].unique())
        ASD = ASD[["gene","Gene stable ID"]] 
        DF = ASD.copy()
    elif gene_lists == "NDD":
        # reading Final NDD genes and make a list out of it
        NDD = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/Gene_lists/genetrek-data-v1-2021-05-18.tsv", sep="\t")
        NDD = NDD[NDD["highConfidenceNddV1"]==True].reset_index(drop=True)
        NDD = NDD.rename(columns={"ensemblID":"Gene stable ID","symbol":"gene"})
        NDD = NDD[["gene","Gene stable ID"]]
        NDD_list = NDD["gene"].to_list()
        Ensemble_list = NDD["Gene stable ID"].to_list()
        # change the name of NDD list to make it compatible iwth the rest of the code
        ASD_list = NDD_list.copy()
        DF = NDD.copy()
    elif gene_lists == "CP":
        # reading Final NDD genes and make a list out of it
        CP = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/Gene_lists/CP_list_updated.tsv", sep="\t")
        CP_list = list(CP["gene"].dropna().unique())
        Ensemble_list = list(CP["Gene stable ID"].dropna().unique())
        CP = CP[["gene","Gene stable ID"]] 
        # change the name of NDD list to make it compatible with the rest of the code
        ASD_list = CP_list.copy()
        DF = CP.copy()
    elif gene_lists == None:
        ASD_list = []
        Ensemble_list = []
        DF = Pd.DataFrame(columns=["gene","Gene stable ID"])
    else:
        print ("please insert the correct gene_lists varibale:NDD or ASD or CP")
        ASD_list = []
        Ensemble_list = []
        DF = Pd.DataFrame(columns=["gene","Gene stable ID"])
        
    return ASD_list, Ensemble_list, DF