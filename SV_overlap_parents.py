#!/usr/bin/env python
#PBS -l walltime=23:59:00
#PBS -l nodes=1:ppn=32
#PBS -l mem=64g,vmem=64g
#PBS -m abe

print ("0")

import pandas as Pd
import numpy as np
from itertools import zip_longest
import dask.dataframe as dd

"""
THIS CODE finds the overlap between CNVs and lncRNAs from lncipedia for all inheritance
PICK the dataset and inheritance and RUN!
"""
print ("1")

DATASET = ["MSSNG","SSC","1000G","MGRB"]
#INHERITANCE = ["P_denovo","Maternal","Paternal"]
CASE = [["proband","affected sibling"],["unaffected sibling","other sibling"],["affected parents"],["unaffaceted parents"]]
Parents = ["father", "mother", "grandfather", "grandmother"]

#pick from above (set case to 0 if you want probands or 1 if you want unaffected)
dataset = DATASET[0]
#inheritance = INHERITANCE[0]
case = CASE[0]

# Some measurements if we pick 1000G or MGRB
if dataset=="1000G":
#    INHERITANCE = [float("NaN")]
    case = "-"
elif dataset=="MGRB":
    case = "singleton"
    


INcRNA = Pd.read_csv("lncipedia_5_2_hg38.bed", sep='\t', header=None)
INcRNA.columns= ["CHR", "Start", "End", "Gene_name", "Score", "Strand", "Start_coding", "end_coding", "Color", "blockCount", "blockSizes", "blockStarts"]
INcRNA["Start"] = INcRNA["Start"]+1

CNV = dd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/SVs."+dataset+".freq_1percent.HQR.tsv", sep='\t',dtype={"Inheritance":"object"})
CNV = CNV[["sample","chrm","start","end","sv_type","Inheritance", "cds_symbol", "gnomAD_pLI","Relation","Affection"]]
CNV_1 = CNV.compute()
CNV_1 = CNV_1.reset_index(drop=True)

print ("2")

#CNV_1 = CNV[CNV.Inheritance.isin(INHERITANCE)]
CNV_l = CNV_1.astype({"start": int, "end": int})

CNV_new = CNV_l.sort_values(['chrm', 'start'], ascending = [1,1]).reset_index(drop=True)

MSSNG = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/"+dataset+"_metadata.tsv", sep='\t')
CNV_new["Relation_updated"] = CNV_new["Relation"]
CNV_new.loc[CNV_new["Relation"]=="child","Relation_updated"] = np.where(CNV_new[CNV_new["Relation"]=="child"]["Affection"]==1,"unaffected sibling","proband")
CNV_new.loc[CNV_new["Relation"].isin(Parents),"Relation_updated"] = np.where(CNV_new[CNV_new["Relation"].isin(Parents)]["Affection"]==1,"unaffaceted parents","affected parents")
#Proband_check = dict(zip(MSSNG["Sample ID"],MSSNG["Relation_updated"]))
#Affection_check = dict(zip(MSSNG["Sample ID"],MSSNG["Affection"]))
Sex_check = dict(zip(MSSNG["Sample ID"],MSSNG["Sex"]))

lncRNA_new = INcRNA.sort_values(['CHR', 'Start'], ascending = [1,1]).reset_index(drop=True)
lncRNA_new["blockSizes"]=lncRNA_new["blockSizes"].str.rstrip(",")
lncRNA_new["blockStarts"]=lncRNA_new["blockStarts"].str.rstrip(",")

lncRNA_new["blocksize_list"] = [list(map(int, i.split(","))) for i in lncRNA_new["blockSizes"]]
lncRNA_new["blockstart_list"] = [list(map(int, i.split(","))) for i in lncRNA_new["blockStarts"]]
lncRNA_new = lncRNA_new.assign(blockend_list=[
    [x + y for x, y in zip_longest(d, c, fillvalue=0)]
    for d, c in zip(lncRNA_new.blockstart_list, lncRNA_new.blocksize_list)
])

chrom = []
for i in range(1,23):
    chrom.append ("chr"+str(i))
chrom.append("chrX")
chrom.append("chrY")

CNV_new_p = CNV_new[CNV_new["Relation_updated"].isin(list(case))].reset_index(drop=True)

A, B, C, D, E, F, G, H, L, M, N  = [],[],[],[],[],[],[],[],[],[],[]
for l in range(len(chrom)):
    lncRNA_new_chr = lncRNA_new[lncRNA_new.CHR == chrom[l]].reset_index(drop=True)
    CNV_chr = CNV_new_p[CNV_new_p.chrm == chrom[l]].reset_index(drop=True)
    for j in range(len(lncRNA_new_chr)):
        CNV_new_chr = CNV_chr[(CNV_chr.start <= lncRNA_new_chr.End[j]) & (CNV_chr.end >= lncRNA_new_chr.Start[j])].reset_index(drop=True)
        for i in range (len(CNV_new_chr)):
            #if Proband_check.get(CNV_new_chr["sample"][i]) in case:
                
            if len(range(max(CNV_new_chr["start"][i], lncRNA_new_chr["Start"][j]), min(CNV_new_chr["end"][i], lncRNA_new_chr["End"][j])+1))>0:

                if ((CNV_new_chr["start"][i] <= lncRNA_new_chr["Start"][j]) and (CNV_new_chr["end"][i] >= lncRNA_new_chr["End"][j])):
                    A.append (CNV_new_chr["sample"][i])
                    B.append (lncRNA_new_chr["Gene_name"][j])
                    C.append("Complete overlap")
                    D.append(chrom[l])
                    E.append(CNV_new_chr["start"][i])
                    F.append(CNV_new_chr["end"][i])
                    G.append(CNV_new_chr["sv_type"][i])
                    H.append(CNV_new_chr["Inheritance"][i])
                    L.append(Sex_check.get(CNV_new_chr["sample"][i]))
                    M.append(CNV_new_chr["cds_symbol"][i])
                    N.append(CNV_new_chr["gnomAD_pLI"][i])
                    print (CNV_new_chr["sample"][i],lncRNA_new_chr["Gene_name"][j],lncRNA_new_chr["CHR"][j])             
                    
                elif ((CNV_new_chr["start"][i] <= lncRNA_new_chr["Start"][j]) and (CNV_new_chr["end"][i] < lncRNA_new_chr["End"][j])):
                    if ((CNV_new_chr["start"][i] <= (lncRNA_new_chr["blockstart_list"][j][0]+lncRNA_new_chr["Start"][j])) and (CNV_new_chr["end"][i] >= (lncRNA_new_chr["blockend_list"][j][0]+lncRNA_new_chr["Start"][j]))):
                        A.append (CNV_new_chr["sample"][i])
                        B.append (lncRNA_new_chr["Gene_name"][j])
                        C.append("Partially overlap")
                        D.append(chrom[l])
                        E.append(CNV_new_chr["start"][i])
                        F.append(CNV_new_chr["end"][i])
                        G.append(CNV_new_chr["sv_type"][i])
                        H.append(CNV_new_chr["Inheritance"][i])
                        L.append(Sex_check.get(CNV_new_chr["sample"][i]))
                        M.append(CNV_new_chr["cds_symbol"][i])
                        N.append(CNV_new_chr["gnomAD_pLI"][i])

                
                
                elif ((CNV_new_chr["start"][i] > lncRNA_new_chr["Start"][j]) and (CNV_new_chr["end"][i] >= lncRNA_new_chr["End"][j])):
                    if ((CNV_new_chr["start"][i] <= (lncRNA_new_chr["blockstart_list"][j][-1]+lncRNA_new_chr["Start"][j])) and (CNV_new_chr["end"][i] >= (lncRNA_new_chr["blockend_list"][j][-1]+lncRNA_new_chr["Start"][j]))):
                        A.append (CNV_new_chr["sample"][i])
                        B.append (lncRNA_new_chr["Gene_name"][j])
                        C.append("Partially overlap")
                        D.append(chrom[l])
                        E.append(CNV_new_chr["start"][i])
                        F.append(CNV_new_chr["end"][i])
                        G.append(CNV_new_chr["sv_type"][i])
                        H.append(CNV_new_chr["Inheritance"][i])
                        L.append(Sex_check.get(CNV_new_chr["sample"][i]))
                        M.append(CNV_new_chr["cds_symbol"][i])
                        N.append(CNV_new_chr["gnomAD_pLI"][i])

                            
                else:
                    for k in range(len(lncRNA_new_chr["blockstart_list"][j])):
                        if ((CNV_new_chr["start"][i] <= (lncRNA_new_chr["blockstart_list"][j][k]+lncRNA_new_chr["Start"][j])) and (CNV_new_chr["end"][i] >= (lncRNA_new_chr["blockend_list"][j][k]+lncRNA_new_chr["Start"][j]))):
                            A.append (CNV_new_chr["sample"][i])
                            B.append (lncRNA_new_chr["Gene_name"][j])
                            C.append("Partially overlap")
                            D.append(chrom[l])
                            E.append(CNV_new_chr["start"][i])
                            F.append(CNV_new_chr["end"][i])
                            G.append(CNV_new_chr["sv_type"][i])
                            H.append(CNV_new_chr["Inheritance"][i])
                            L.append(Sex_check.get(CNV_new_chr["sample"][i]))
                            M.append(CNV_new_chr["cds_symbol"][i])
                            N.append(CNV_new_chr["gnomAD_pLI"][i])


                            print (k,CNV_new_chr["sample"][i],lncRNA_new_chr["Gene_name"][j],lncRNA_new_chr["CHR"][j])
                            break

                        else:
                            A.append (CNV_new_chr["sample"][i])
                            B.append (lncRNA_new_chr["Gene_name"][j])
                            C.append("Partially exon")
                            D.append(chrom[l])
                            E.append(CNV_new_chr["start"][i])
                            F.append(CNV_new_chr["end"][i])
                            G.append(CNV_new_chr["sv_type"][i])
                            H.append(CNV_new_chr["Inheritance"][i])
                            L.append(Sex_check.get(CNV_new_chr["sample"][i]))
                            M.append(CNV_new_chr["cds_symbol"][i])
                            N.append(CNV_new_chr["gnomAD_pLI"][i])

                               
                
            
                                                
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



if case==CASE[0]:
    caseout = ""
elif case == CASE[1]:
    caseout = "unaffected"
elif case == CASE[2]:
    caseout = "parents"
else:
    caseout = ""
df.to_csv("/hpf/largeprojects/tcagstor/users/sghaffari/lncRNA_Overlaps_IO/lncRNA_Overlapped_SV_"+dataset+"_"+caseout+".tsv", index=False, sep='\t')                    
    

