#!/usr/bin/env python
#PBS -l walltime=23:59:00
#PBS -l nodes=1:ppn=32
#PBS -l mem=64g,vmem=64g
#PBS -m abe

import pandas as Pd
import numpy as np
from itertools import zip_longest

"""
THIS CODE finds the overlap between Copy Number Variations (CNVs) and lncRNAs from lncipedia database for all inheritance
PICK the dataset and inheritance and RUN!
"""

DATASET = ["MSSNG","SSC","1000G","MGRB"]
#INHERITANCE = ["P_denovo","Maternal","Paternal"]
CASE = [["proband","affected sibling"],["unaffected sibling","other sibling"],["affected parents"],["unaffaceted fathers"],["unaffaceted parents"]]
Parents = ["father", "mother", "grandfather", "grandmother"]


# pick manually from above (set case to 0 if you want probands or 1 if you want unaffected, 2 for affected parents, 3 unaffacted fathers, 4 unaffacted parents)
dataset = DATASET[0]
# inheritance = INHERITANCE[0]
case = CASE[2]

# load LNCIPEDIA database hg38 bed file for lncRNAs
INcRNA = Pd.read_csv("/home/sghaffari/lncipedia_5_2_hg38.bed", sep='\t', header=None)
INcRNA.columns= ["CHR", "Start", "End", "Gene_name", "Score", "Strand", "Start_coding", "end_coding", "Color", "blockCount", "blockSizes", "blockStarts"]
INcRNA["Start"] = INcRNA["Start"]+1

lncRNA_new = INcRNA.sort_values(['CHR', 'Start'], ascending = [1,1]).reset_index(drop=True)
lncRNA_new["blockSizes"]=lncRNA_new["blockSizes"].str.rstrip(",")
lncRNA_new["blockStarts"]=lncRNA_new["blockStarts"].str.rstrip(",")

lncRNA_new["blocksize_list"] = [list(map(int, i.split(","))) for i in lncRNA_new["blockSizes"]]
lncRNA_new["blockstart_list"] = [list(map(int, i.split(","))) for i in lncRNA_new["blockStarts"]]
lncRNA_new = lncRNA_new.assign(blockend_list=[
            [x + y for x, y in zip_longest(d, c, fillvalue=0)]
            for d, c in zip(lncRNA_new.blockstart_list, lncRNA_new.blocksize_list)
            ])


# Automatically iterate over all neede cases and datasets 
for data_i in range(4):
    dataset = DATASET[data_i]
    if (data_i<2): CASE_NUM = 4
    else : CASE_NUM = 1
    for case_j in range(CASE_NUM):
        case = CASE[case_j]

        # Some measurements if we pick 1000G or MGRB
        if dataset=="1000G":
            #    INHERITANCE = [float("NaN")]
            case = ["-"]
        elif dataset=="MGRB":
            case = ["singleton"]
    
        
        # load CNV database - make sure it is in hg38 and the column names should be in a format of "START", "END" , "CHROM" for beginning, end and the chromosome of the CNVs loci
        CNV = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/CNVs."+dataset+".freq_1percent.HQR.tsv", sep='\t',low_memory=False)
        #CNV_1 = CNV[CNV.Inheritance.isin(INHERITANCE)]
        CNV_l = CNV.astype({"START": int, "END": int})
        
        CNV_new = CNV_l.sort_values(['CHROM', 'START'], ascending = [1,1]).reset_index(drop=True)
        
        # load the database of patients information
        MSSNG = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/"+dataset+"_metadata.tsv", sep='\t')
        MSSNG["Relation_updated"] = MSSNG["Relation"]
        MSSNG.loc[MSSNG["Relation"]=="child","Relation_updated"] = np.where(MSSNG[MSSNG["Relation"]=="child"]["Affection"]==1,"unaffected sibling","proband")
        MSSNG.loc[MSSNG["Relation"].isin(Parents),"Relation_updated"] = np.where(MSSNG[MSSNG["Relation"].isin(Parents)]["Affection"]==1,"unaffaceted parents","affected parents")
        MSSNG.loc[MSSNG["Relation"]=="father","Relation_updated"] = np.where(MSSNG[MSSNG["Relation"]=="father"]["Affection"]==1,"unaffaceted fathers","affected parents")
        Proband_check = dict(zip(MSSNG["Sample ID"],MSSNG["Relation_updated"]))
        #Affection_check = dict(zip(MSSNG["Sample ID"],MSSNG["Affection"]))
        #Sex_check = dict(zip(MSSNG["Sample ID"],MSSNG["Sex"]))
        
        # merge Relation column to our CNV database
        MSSNG_updated = MSSNG[["Sample ID","Relation_updated","Sex","Predicted ancestry"]]
        MSSNG_updated = MSSNG_updated.rename(columns={"Sample ID":"sample"})
        CNV_new2 = CNV_new.merge(MSSNG_updated, 
                         on="sample",
                         how="left")
    
    
        # We care about unffacted fathers only in X-Chromosome
        if case == CASE[3]:
            CNV_new2 = CNV_new2[CNV_new2["CHROM"]=="chrX"]

        # Filter a database that is only in our desired case
        CNV_new_p = CNV_new2[CNV_new2["Relation_updated"].isin(case)].reset_index(drop=True)
        
        chrom = []
        for i in range(1,23):
            chrom.append ("chr"+str(i))
        chrom.append("chrX")
        chrom.append("chrY")

        A, B, C, D, E, F, G, H, L, M, N, O  = [],[],[],[],[],[],[],[],[],[],[],[]
        for l in range(len(chrom)):
            lncRNA_new_chr = lncRNA_new[lncRNA_new.CHR == chrom[l]].reset_index(drop=True)
            CNV_chr = CNV_new_p[CNV_new_p.CHROM == chrom[l]].reset_index(drop=True)
            for j in range(len(lncRNA_new_chr)):
                CNV_new_chr = CNV_chr[(CNV_chr.START <= lncRNA_new_chr.End[j]) & (CNV_chr.END >= lncRNA_new_chr.Start[j])].reset_index(drop=True)
                for i in range (len(CNV_new_chr)):
                    if Proband_check.get(CNV_new_chr["sample"][i]) in case:
                    
                        if len(range(max(CNV_new_chr["START"][i], lncRNA_new_chr["Start"][j]), min(CNV_new_chr["END"][i], lncRNA_new_chr["End"][j])+1))>0:
                            if ((CNV_new_chr["START"][i] <= lncRNA_new_chr["Start"][j]) and (CNV_new_chr["END"][i] >= lncRNA_new_chr["End"][j])):
                                A.append (CNV_new_chr["sample"][i])
                                B.append (lncRNA_new_chr["Gene_name"][j])
                                C.append("Complete overlap")
                                D.append(chrom[l])
                                E.append(CNV_new_chr["START"][i])
                                F.append(CNV_new_chr["END"][i])
                                G.append(CNV_new_chr["SVTYPE"][i])
                                H.append(CNV_new_chr["Inheritance"][i])
                                L.append(CNV_new_chr["Sex"][i])
                                M.append(CNV_new_chr["cds_symbol"][i])
                                N.append(CNV_new_chr["gnomAD_pLI"][i])
                                O.append(CNV_new_chr["Predicted ancestry"][i]) 
                            
                            elif ((CNV_new_chr["START"][i] <= lncRNA_new_chr["Start"][j]) and (CNV_new_chr["END"][i] < lncRNA_new_chr["End"][j])):
                                if ((CNV_new_chr["START"][i] <= (lncRNA_new_chr["blockstart_list"][j][0]+lncRNA_new_chr["Start"][j])) and (CNV_new_chr["END"][i] >= (lncRNA_new_chr["blockend_list"][j][0]+lncRNA_new_chr["Start"][j]))):
                                    A.append (CNV_new_chr["sample"][i])
                                    B.append (lncRNA_new_chr["Gene_name"][j])
                                    C.append("Partially overlap")
                                    D.append(chrom[l])
                                    E.append(CNV_new_chr["START"][i])
                                    F.append(CNV_new_chr["END"][i])
                                    G.append(CNV_new_chr["SVTYPE"][i])
                                    H.append(CNV_new_chr["Inheritance"][i])
                                    L.append(CNV_new_chr["Sex"][i])
                                    M.append(CNV_new_chr["cds_symbol"][i])
                                    N.append(CNV_new_chr["gnomAD_pLI"][i])
                                    O.append(CNV_new_chr["Predicted ancestry"][i])
        
                        
                        
                            elif ((CNV_new_chr["START"][i] > lncRNA_new_chr["Start"][j]) and (CNV_new_chr["END"][i] >= lncRNA_new_chr["End"][j])):
                                if ((CNV_new_chr["START"][i] <= (lncRNA_new_chr["blockstart_list"][j][-1]+lncRNA_new_chr["Start"][j])) and (CNV_new_chr["END"][i] >= (lncRNA_new_chr["blockend_list"][j][-1]+lncRNA_new_chr["Start"][j]))):
                                    A.append (CNV_new_chr["sample"][i])
                                    B.append (lncRNA_new_chr["Gene_name"][j])
                                    C.append("Partially overlap")
                                    D.append(chrom[l])
                                    E.append(CNV_new_chr["START"][i])
                                    F.append(CNV_new_chr["END"][i])
                                    G.append(CNV_new_chr["SVTYPE"][i])
                                    H.append(CNV_new_chr["Inheritance"][i])
                                    L.append(CNV_new_chr["Sex"][i])
                                    M.append(CNV_new_chr["cds_symbol"][i])
                                    N.append(CNV_new_chr["gnomAD_pLI"][i])
                                    O.append(CNV_new_chr["Predicted ancestry"][i])
        
                                    
                            else:
                                    for k in range(len(lncRNA_new_chr["blockstart_list"][j])):
                                        if ((CNV_new_chr["START"][i] <= (lncRNA_new_chr["blockstart_list"][j][k]+lncRNA_new_chr["Start"][j])) and (CNV_new_chr["END"][i] >= (lncRNA_new_chr["blockend_list"][j][k]+lncRNA_new_chr["Start"][j]))):
                                            A.append (CNV_new_chr["sample"][i])
                                            B.append (lncRNA_new_chr["Gene_name"][j])
                                            C.append("Partially overlap")
                                            D.append(chrom[l])
                                            E.append(CNV_new_chr["START"][i])
                                            F.append(CNV_new_chr["END"][i])
                                            G.append(CNV_new_chr["SVTYPE"][i])
                                            H.append(CNV_new_chr["Inheritance"][i])
                                            L.append(CNV_new_chr["Sex"][i])
                                            M.append(CNV_new_chr["cds_symbol"][i])
                                            N.append(CNV_new_chr["gnomAD_pLI"][i])
                                            O.append(CNV_new_chr["Predicted ancestry"][i])
        
        
                                            print (k,CNV_new_chr["sample"][i],lncRNA_new_chr["Gene_name"][j],lncRNA_new_chr["CHR"][j])
                                            break                               
                        
                        
                    
        df = Pd.DataFrame()
        df ["Sample"] = A
        df ["lncRNA"] = B
        df ["Description"] = C
        df ["CHROM"] = D
        df ["Start"] = E
        df ["End"] = F
        df ["SVTYPE"] = G
        df ["Inheritance"] = H
        df ["SEX"] = L
        df ["cds_symbol"] = M
        df ["gnomAD_pLI"] = N
        df ["Predicted_ancestry"] = O
        



        if case==CASE[0]:
            caseout = ""
        elif case == CASE[1]:
            caseout = "unaffected"
        elif case == CASE[2]:
            caseout = "parents"
        elif case == CASE[3]:
            caseout = "fathers"
        else:
            caseout = ""
        df.to_csv("/hpf/largeprojects/tcagstor/users/sghaffari/lncRNA_Overlaps_IO/lncRNA_Overlapped_"+dataset+"_"+caseout+".tsv", index=False, sep='\t')                    
    

