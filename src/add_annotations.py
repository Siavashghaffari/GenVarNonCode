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


def ADD_annotations (df, Overlap_probands, Overlap_unaffected, Overlap_parents, case):

    """
    this function get the final databse and add a column of lncRNA annotation to our original MSSNG and SSC databse
    """
    lncRNA_list = df["lncipediaGeneID"].tolist()
    CNV_MSSNG = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/CNVs.MSSNG.freq_1percent.HQR.tsv", sep='\t', low_memory=False)
    MSSNG = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/MSSNG_metadata.tsv", sep='\t')
    CNV_SSC = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/CNVs.SSC.freq_1percent.HQR.tsv", sep='\t', low_memory=False)
    SSC = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/SSC_metadata.tsv", sep='\t')
    
    MSSNG = MSSNG.rename(columns={"Sample ID":"sample"})
    SSC = SSC.rename(columns={"Sample ID":"sample"})
    
    MSSNG_2 = MSSNG[["sample","Sex","Relation","Affection","Mother ID","Father ID","Family ID","Predicted ancestry","Family size"]]
    SSC_2 = SSC[["sample","Sex","Relation","Affection","Mother ID","Father ID","Family ID","Predicted ancestry","Family size"]]

    MSSNG_final = (
    CNV_MSSNG.merge(MSSNG_2, 
              on=['sample'],
              how='left'))
    
    SSC_final = (
    CNV_SSC.merge(SSC_2, 
              on=['sample'],
              how='left'))
    
    CNV = Pd.concat([MSSNG_final,SSC_final])

    # Annot_p is a subset of overlap_probands to be meeged with CNV file and Annot_u is a a subset of overlap_unaffected to be merged with CNV file
    Annot_p = Overlap_probands[Overlap_probands.lncRNA_gene.isin(lncRNA_list)].reset_index(drop=True)
    Annot_u = Overlap_unaffected[Overlap_unaffected.lncRNA_gene.isin(lncRNA_list)].reset_index(drop=True)
    Annot_pa = Overlap_parents[Overlap_parents.lncRNA_gene.isin(lncRNA_list)].reset_index(drop=True)
    
    Annot_p = Annot_p[["Sample","lncRNA_gene","CHROM","Start","End","SVTYPE"]]
    Annot_u = Annot_u[["Sample","lncRNA_gene","CHROM","Start","End","SVTYPE"]]
    Annot_pa = Annot_pa[["Sample","lncRNA_gene","CHROM","Start","End","SVTYPE"]]
    
    Annot_p = Annot_p.drop_duplicates().reset_index(drop=True)
    Annot_u = Annot_u.drop_duplicates().reset_index(drop=True)
    Annot_pa = Annot_pa.drop_duplicates().reset_index(drop=True)
    
    Annot_p = Annot_p.rename(columns={"Sample":"sample","Start":"START","End":'END'})
    Annot_u = Annot_u.rename(columns={"Sample":"sample","Start":"START","End":'END'})
    Annot_pa = Annot_pa.rename(columns={"Sample":"sample","Start":"START","End":'END'})
    
    
    CNV_2 = CNV[["sample","CHROM","START","END","SVTYPE"]]
    
    # dff_p, dfff_p and Final_p are related to probands, and dff_u, dfff_u and Final_u are related to unaffected, dff_pa, dfff_pa and Final_pa are related to affected parents
    dff_p = (
    CNV_2.merge(Annot_p, 
              on=['sample',"CHROM","START","END","SVTYPE"],
              how='right'))
    dff_p = dff_p[dff_p["lncRNA_gene"].notna()]
    dfff_p = dff_p.groupby(["sample","CHROM","START","END","SVTYPE"]).lncRNA_gene.agg(list)
    dfff_p = dfff_p.reset_index(drop=False)
    dfff_p = dfff_p.rename(columns={"lncRNA_gene":"lncipediaGeneID"})
    
    Final_p = (
    CNV.merge(dfff_p, 
              on=['sample','CHROM','START','END','SVTYPE'],
              how='right'))
    
    dff_u = (
    CNV_2.merge(Annot_u, 
              on=['sample',"CHROM","START","END","SVTYPE"],
              how='right'))
    dff_u = dff_u[dff_u["lncRNA_gene"].notna()]
    dfff_u = dff_u.groupby(["sample","CHROM","START","END","SVTYPE"]).lncRNA_gene.agg(list)
    dfff_u = dfff_u.reset_index(drop=False)
    dfff_u = dfff_u.rename(columns={"lncRNA_gene":"lncipediaGeneID"})
    Final_u = (
    CNV.merge(dfff_u, 
              on=['sample','CHROM','START','END','SVTYPE'],
              how='right'))

    dff_pa = (
    CNV_2.merge(Annot_pa, 
              on=['sample',"CHROM","START","END","SVTYPE"],
              how='right'))
    dff_pa = dff_pa[dff_pa["lncRNA_gene"].notna()]
    dfff_pa = dff_pa.groupby(["sample","CHROM","START","END","SVTYPE"]).lncRNA_gene.agg(list)
    dfff_pa = dfff_pa.reset_index(drop=False)
    dfff_pa = dfff_pa.rename(columns={"lncRNA_gene":"lncipediaGeneID"})  
    Final_pa = (
    CNV.merge(dfff_pa, 
              on=['sample','CHROM','START','END','SVTYPE'],
              how='right'))
    
    # Combine probands and unaffected data
    DF_FINAL = Pd.concat([Final_p, Final_u, Final_pa])
    # Publish the final CNV database with added annotation
    DF_FINAL.to_csv("Final_ASD_overlaps_info"+case+".tsv", sep="\t")
    
    
    # Check for having multiple samples from the same family
    FAM_p = (
    CNV.merge(Annot_p, 
              on=['sample',"CHROM","START","END","SVTYPE"],
              how='right'))
    FAM_p = FAM_p[dff_p["lncRNA_gene"].notna()]
    FAM_p = FAM_p.reset_index(drop=True)
    FAM_p = FAM_p.rename(columns={"lncRNA_gene":"lncipediaGeneID"})
    
    FAM_u = (
    CNV.merge(Annot_u, 
              on=['sample',"CHROM","START","END","SVTYPE"],
              how='right'))
    FAM_u = FAM_u[dff_u["lncRNA_gene"].notna()]
    FAM_u = FAM_u.reset_index(drop=True)
    FAM_u = FAM_u.rename(columns={"lncRNA_gene":"lncipediaGeneID"})
    
    FAM_pa = (
    CNV.merge(Annot_pa, 
              on=['sample',"CHROM","START","END","SVTYPE"],
              how='right'))
    FAM_pa = FAM_pa[dff_pa["lncRNA_gene"].notna()]
    FAM_pa = FAM_pa.reset_index(drop=True)
    FAM_pa = FAM_pa.rename(columns={"lncRNA_gene":"lncipediaGeneID"})
    
    # Combine probands and unaffected data
    FAM_FINAL = Pd.concat([FAM_p, FAM_u, FAM_pa])
    FAM_2 = FAM_FINAL.groupby(["lncipediaGeneID","SVTYPE"])["Family ID"].agg(list).reset_index(drop=False)
    FAM_2["Same_Family"] = FAM_2["Family ID"].apply(lambda x:list(set([a for a in x if x.count(a)>1])))
    FAM_2["N_Family"] = FAM_2["Same_Family"].apply(lambda x:len(x))
    FAM_2 = FAM_2[["lncipediaGeneID","SVTYPE","Same_Family","N_Family"]]
    
    # Keep only families with multiple samples
    FAM_3 = FAM_2[FAM_2.N_Family>0].reset_index(drop = True)
    FAM_4 = FAM_3.explode('Same_Family', ignore_index = True)
    FAM_5 = FAM_4.drop_duplicates(["lncipediaGeneID","SVTYPE","Same_Family"]).reset_index(drop=True)
    FAM_6 = list(zip(FAM_5["lncipediaGeneID"], FAM_5["SVTYPE"], FAM_5["Same_Family"]))
    
    # Binarize affection and Sex 
    FAM_FINAL["Affection_binary"] = np.where(FAM_FINAL["Affection"].isin([0,2]),1,0)
    #FAM_FINAL["Affection_binary"] = np.where(FAM_FINAL["Affection"]==1,0,1)
    FAM_FINAL["Sex_binary"] = np.where(FAM_FINAL["Sex"]=="male",1,0)
    #FAM_FINAL["Sex_binary"] = np.where(FAM_FINAL["Sex"]=="female",0,1)
    
    
    # Organizing data from the same family
    A, B, C, D, E, F = [], [], [], [], [], []
    for a,b,c in FAM_6:
        D.append(FAM_FINAL[(FAM_FINAL["lncipediaGeneID"]==a) & (FAM_FINAL["SVTYPE"]==b) & (FAM_FINAL["Family ID"]==c)]["Relation"].to_list())
        E.append(FAM_FINAL[(FAM_FINAL["lncipediaGeneID"]==a) & (FAM_FINAL["SVTYPE"]==b) & (FAM_FINAL["Family ID"]==c)]["Affection_binary"].to_list())
        F.append(FAM_FINAL[(FAM_FINAL["lncipediaGeneID"]==a) & (FAM_FINAL["SVTYPE"]==b) & (FAM_FINAL["Family ID"]==c)]["Sex_binary"].to_list())
        A.append(a)
        B.append(b)
        C.append(c)
    
    FAM_7 = Pd.DataFrame()
    FAM_7 ["lncipediaGeneID"] = A
    FAM_7 ["SVTYPE"] = B
    FAM_7 ["Same_Family"] = C
    FAM_7 ["Relation_F"] = D
    FAM_7 ["Affection_F"] = E
    FAM_7 ["Sex_F"] = F
    
    FAM_8 = FAM_7.groupby(["lncipediaGeneID","SVTYPE"]).agg({"Same_Family": lambda x: x.tolist(),"Relation_F":lambda x: x.tolist(),
                                                          "Affection_F": lambda x: x.tolist(),"Sex_F":lambda x: x.tolist()}).reset_index()
    FAM_9 = FAM_8[["lncipediaGeneID","SVTYPE","Relation_F", "Affection_F", "Sex_F"]]
    
    
    
    # Sanity Check for DGV database
    df_Final_test = (
    CNV.merge(dff_p, 
              on=['sample',"START","END","SVTYPE"],
              how='left'))
    
    # DGV_coverage_freq
    DGV1 = df_Final_test.groupby(["lncRNA_gene","SVTYPE"])["DGVpercFreq_subjects_coverageStudies_50percRecipOverlap"].mean()
    DGV_freq1 = DGV1.to_frame().reset_index()
    DGV_freq1 = DGV_freq1.rename(columns={"lncRNA_gene":"lncipediaGeneID","DGVpercFreq_subjects_coverageStudies_50percRecipOverlap":"DGV Freq"})
    # DGV_number_studies
    DGV2 = df_Final_test.groupby(["lncRNA_gene","SVTYPE"])["DGV_N_studies_50percRecipOverlap"].mean()
    DGV_freq2 = DGV2.to_frame().reset_index()
    DGV_freq2 = DGV_freq2.rename(columns={"lncRNA_gene":"lncipediaGeneID","DGV_N_studies_50percRecipOverlap":"DGV N_Studies"})
    # Merge above databases
    DGV_freq = (
    DGV_freq1.merge(DGV_freq2, 
              on=["lncipediaGeneID","SVTYPE"],
              how='left'))
        
    df_DGV_test = (
    df.merge(DGV_freq, 
              on=["lncipediaGeneID","SVTYPE"],
              how='left'))
    
    # Add family info to our final DB 
    df_DGV_test = (
    df_DGV_test.merge(FAM_2, 
              on=["lncipediaGeneID","SVTYPE"],
              how='left'))
    
    df_DGV_test = (
    df_DGV_test.merge(FAM_9, 
              on=["lncipediaGeneID","SVTYPE"],
              how='left'))
    
    #df_DGV_test = df_DGV_test[["lncipediaGeneID","Exon_no","SVTYPE","DGV Freq"]]
    df_DGV_test["DGV Freq"] = df_DGV_test["DGV Freq"].fillna(0)
    df_DGV_test["DGV N_Studies"] = df_DGV_test["DGV N_Studies"].fillna(0)
    
    return df_DGV_test
