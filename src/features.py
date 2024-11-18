#!/usr/bin/env python
#PBS -l walltime=47:59:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=64g,vmem=64g
#PBS -m e

"""
Created on Tue Aug 31 10:42:26 2021

@author: Siavash Ghaffari
"""

import pandas as Pd
import numpy as np
from gene_list_pick import Gene_list_pick


class Features(object):
    """
    Features class encapsulates all the logic necessary for evaluating the importance of the lncRNA_list.
    The Features accepts lncRNA exons overlapped with CNVs data and find the ones that are belived to have some kind of importance.
    You will need to first construct a Features instance, passing the dataset to the constructor. You will then call the methods to 
    annotate different fetaures on the overlapped lncRNAs. 
    """
    def __init__(self, Overlap, **kwargs):
        """
        Construct a new Feature instance.
        Required arguments:
        - Overlap: a list of final overalpped lncRNAs exons by CNVs
        - feat: one of features 
        """
        
        self.Overlap = Overlap
        #GERP = Pd.read_csv("C:/Siavash/Conservation Data/Conserved_elements_hg38.tsv", sep='\t')
        #VISTA = Pd.read_csv("C:/Siavash/Conservation Data/VISTA_Enhancer_hg38.tsv", sep='\t')
        #Ultra = Pd.read_csv("C:/Siavash/Conservation Data/Ultra_conserved_hg38.tsv", sep='\t')

        
        #self.Overlap = Pd.read_csv("Exons_probands+unaffected.tsv", sep='\t',low_memory=False,converters={"Exon_Co": eval,"Exons": eval})
    
    def Conserved_exons (self, Cons, if_Score=False, if_Enhancer=False, if_intercation=False):
        """
        Find overlaps between Exons covered by CNVs and Conserved elemnets
        Cons is conserved_elements dataframe coming from different databases
        if _Score is a boolean variable indicates if to create the score column
        if _Enhancer is a boolean variable indicates if to create the enhancer column
        """
        

        self.Overlap["Start_exon"] = self.Overlap["Exon_Co"].map(lambda x:x[0])
        self.Overlap["End_exon"] = self.Overlap["Exon_Co"].map(lambda x:x[-1])

        # Overlap Algorithm
        chrom = []
        for l in range(1,23):
            chrom.append ("chr"+str(l))
        chrom.append("chrX")
        chrom.append("chrY")

        A, B, C, D, E, F, G, H, I, J, K, L, M, N  = [],[],[],[],[],[],[],[],[],[],[],[],[],[]
        for l in range(len(chrom)):
            Cons_chr = Cons[Cons.Start_chr == chrom[l]].reset_index(drop=True)
            Overlap_chr = self.Overlap[self.Overlap["CHROM"] == chrom[l]].reset_index(drop=True)
            for i in range(len(Cons_chr)):
                Overlap_chr_range =  Overlap_chr[(Overlap_chr.Start_exon > Cons_chr.Start_hg38[i]-2000000) & (Overlap_chr.End_exon < Cons_chr.End_hg38[i]+2000000)].reset_index(drop=True)
                for j in range(len(Overlap_chr_range)):
                    inter = len(range(max(Cons_chr["Start_hg38"][i], Overlap_chr_range["Start_exon"][j]), min(Cons_chr["End_hg38"][i], Overlap_chr_range["End_exon"][j])))
                    if inter>0:
                        A.append (Overlap_chr_range["Exon_no"][j])
                        B.append (Overlap_chr_range["lncipediaGeneID"][j])
                        C.append (Overlap_chr_range["CHROM"][j])
                        D.append(Overlap_chr_range["SVTYPE"][j])
                        E.append(Overlap_chr_range["Exon_Co"][j])
                        F.append(inter)
                        G.append((Cons_chr["Start_hg38"][i],Cons_chr["End_hg38"][i]))
                        if if_Enhancer:
                            H.append(Cons_chr.Enhancer[i])
                        if if_Score:
                            H.append(Cons_chr.Score[i])
                        if if_intercation:
                            H.append(Cons_chr.Gene_ID[i])
                            I.append(Cons_chr.Gene_Name[i])
                            J.append(Cons_chr.Gene_CHR[i])
                            K.append(Cons_chr.Gene_TSS[i])
                            L.append(Cons_chr.Gene_Strand[i])
                            M.append(Cons_chr.Score[i]) 
                        
                        #print (Cons_chr["Start_hg38"][i],Cons_chr["End_hg38"][i],Overlap_chr_range["Exon_no"][j] , chrom[l],Overlap_chr_range["lncipediaGeneID"][j],inter) 
                
        Conserved_exons = Pd.DataFrame()
        Conserved_exons ["lncRNA"] = B
        Conserved_exons ["Exon_no"] = A
        Conserved_exons ["CHROM"] = C
        Conserved_exons ["SVTYPE"] = D
        Conserved_exons ["Exons_Co"] = E
        Conserved_exons ["length"] = F
        Conserved_exons ["elements_Co"] = G
        if if_Enhancer:
            Conserved_exons ["Enhancer"] = H
        if if_Score:
            Conserved_exons ["Score"] = H
        if if_intercation:
            Conserved_exons ["Gene_ID"] = H
            Conserved_exons ["Gene_Name"] = I
            Conserved_exons ["Gene_CHR"] = J
            Conserved_exons ["Gene_TSS"] = K
            Conserved_exons ["Gene_Strand"] = L
            Conserved_exons ["Score"] = M    
             
        #Conserved_exons.to_csv("Conserved_Exons_genes.tsv", sep="\t")
        return Conserved_exons
    
    def Make_Conversion_list_for_Coexp (self, LncBook, lncRNAKB, GTEX):
        """
        this function convert lncRNA list from LNCPEDIA to LncBook and lncRNAKB
        """
        lncRNA_list = list(self.Overlap["lncipediaGeneID"].unique())
        # CALL the lncRNA_Convertor method
        DATASET = [LncBook,lncRNAKB]
        LncBook_Convertor = self.lncRNA_Convertor (DATASET[0], self.Overlap, "LncBook")
        lncRNAKB_Convertor = self.lncRNA_Convertor (DATASET[1], self.Overlap, "lncRNAKB")
        
        # finalizing and cleaning LncBook_Convertor
        LncBook_Convertor_2 = LncBook_Convertor.drop_duplicates().reset_index(drop=True)
        LncBook_Convertor_3 = LncBook_Convertor_2.loc[LncBook_Convertor_2.groupby(['Lincpedia','LncBook'])['inter'].idxmax()].reset_index(drop=True)
        Low = LncBook_Convertor_3[~((LncBook_Convertor_3.overlap_1>0.99) | (LncBook_Convertor_3.overlap_2>0.99))]
        Low=Low[Low.overlap<2]
        Low = Low[Low.inter<0.8]
        self.LncBook_Final = Pd.concat([LncBook_Convertor_3, Low, Low]).drop_duplicates(keep=False)
        self.LncBook_final = self.LncBook_Final[["Lincpedia","LncBook"]]
        # make a list LncBook_Final
        self.LncBook_Final_list = list(self.LncBook_Final["LncBook"].unique())
        
        # finalizing and cleaning lncRNAKB_Convertor
        lncRNAKB_Convertor_2 = lncRNAKB_Convertor.drop_duplicates().reset_index(drop=True)
        lncRNAKB_Convertor_3 = lncRNAKB_Convertor_2.loc[lncRNAKB_Convertor_2.groupby(['Lincpedia','lncRNAKB'])['inter'].idxmax()].reset_index(drop=True)   
        Low = lncRNAKB_Convertor_3[~((lncRNAKB_Convertor_3.overlap_1>0.99) | (lncRNAKB_Convertor_3.overlap_2>0.99))]
        Low=Low[Low.overlap<2]
        Low=Low[Low.inter<0.8]
        lncRNAKB_Final = Pd.concat([lncRNAKB_Convertor, Low, Low]).drop_duplicates(keep=False)
        # Import GTEX conversion
        GTEX = GTEX[["lncipediaGeneID","GENE_symbol"]]
        GTEX = GTEX.drop_duplicates().reset_index(drop=True)
        GTEX = GTEX.rename(columns={"lncipediaGeneID":"Lincpedia"})
        lncRNAKB_Final_updated = (
            lncRNAKB_Final.merge(GTEX, 
              on=['Lincpedia'],
              how='left'))
        self.lncRNAKB_Final_updated = lncRNAKB_Final_updated.drop_duplicates().reset_index(drop=True)
        self.lncRNAKB_Final_updated = self.lncRNAKB_Final_updated[["Lincpedia","lncRNAKB","GENE_symbol"]]
        self.LncRNAKB_Final_list_A = list(self.lncRNAKB_Final_updated["lncRNAKB"].unique())
        self.LncRNAKB_Final_list_B = list(self.lncRNAKB_Final_updated["Lincpedia"].unique())
        self.LncRNAKB_Final_list_C = list(self.lncRNAKB_Final_updated["GENE_symbol"].unique())
        
    # Convertor function based on information from overlapped exon    
    def lncRNA_Convertor (self, Database, Exons, name):
        """
        Convertor function based on information from overlapped exon
        """
        lncRNA_list = list(Exons["lncipediaGeneID"].unique())
        A, B, C, D, E, F, G, H, I, J, K, L = [],[],[],[],[],[],[],[],[],[],[],[]
        for i in range(len(lncRNA_list)):
            Exon_lncRNA = Exons[Exons["lncipediaGeneID"]==lncRNA_list[i]].reset_index(drop=True)
            Database_lncRNA = Database[Database["Lincpedia"]==lncRNA_list[i]].reset_index(drop=True)
            for j in range(len(Exon_lncRNA)):
                for k in range(len(Database_lncRNA)):
                    inter = len(range(max(Exon_lncRNA["Exon_Co"][j][0], Database_lncRNA["Start_"+name][k]), min(Exon_lncRNA["Exon_Co"][j][1], Database_lncRNA["End_"+name][k])))
                    if inter>0:
                        A.append(Database_lncRNA["Lincpedia"][k])
                        B.append(Database_lncRNA[name][k])
                        C.append(Database_lncRNA["overlap_1"][k])
                        D.append(Database_lncRNA["overlap_2"][k])
                        E.append(Database_lncRNA["overlap"][k])
                        F.append(inter/(Exon_lncRNA["Exon_Co"][j][1]-Exon_lncRNA["Exon_Co"][j][0]))
                    
        Convertor = Pd.DataFrame()
        Convertor ["Lincpedia"] = A
        Convertor [name] = B
        Convertor ["overlap_1"] = C
        Convertor ["overlap_2"] = D
        Convertor ["overlap"] = E
        Convertor ["inter"] = F
    
        return Convertor
    
                    
                                    
    def Co_exp_LncBook (self, coexp, gene_list):
        """
        this function returns lncRNAs from our lists that are co-expressed with ASD or NDD or CP gene list
        coexp is a database for each cell_line in LncBook
        gene_list is a string "ASD", "NDD", "CP"
        """
        
        coexp.columns= ["Gene1","Gene2","r","p_value","connectivity","Kdiff"]
        coexp["Gene2"] = coexp["Gene2"].str.split(".").map(lambda x:x[0])
        coexp["Gene1"] = coexp["Gene1"].str.split(".").map(lambda x:x[0])

        A = self.LncBook_Final_list.copy() 
        B = Gene_list_pick (gene_list)[1]
        Final = coexp[coexp["Gene1"].isin(A)&coexp["Gene2"].isin(B)]
        Final2 = coexp[coexp["Gene1"].isin(B)&coexp["Gene2"].isin(A)]
        Final3 = Pd.concat([Final,Final2])
        #set the P-value cut-off
        Final4 = Final3[Final3["p_value"]<=0.01]
        #Final4 = Final3.copy()
        Final4["Status"] = Final4.Gene1.str.contains("HSALNG")
        Final4["LncBook"] = [d if e else c for d, c, e in zip(Final4.Gene1,Final4.Gene2,Final4.Status)]
        Final4["Gene stable ID"] = [c if e else d for d, c, e in zip(Final4.Gene1,Final4.Gene2,Final4.Status)]
        #ASD = ASD.rename(columns={"Gene stable ID":"ASD_ensemble", "gene":"ASD Gene"})
        Genes = Gene_list_pick (gene_list)[2]
        Final5 = (
            Final4.merge(Genes, 
            on=['Gene stable ID'],
              how='left'))
        Final6 = (
            Final5.merge(self.LncBook_final, 
              on=['LncBook'],
              how='left'))
        #Final_6.columns
        return Final6
        
                 
    def Co_exp_lncrnakb (self, gene_list):
        """
        this function returns lncRNAs from our lists that are co-expressed with ASD or NDD or CP gene list
        coexp is a database for each cell_line in lncrnakb
        gene_list is a string "ASD", "NDD", "CP"
        """
        
        A = self.LncRNAKB_Final_list_A.copy()
        B = self.LncRNAKB_Final_list_B.copy()
        C = self.LncRNAKB_Final_list_C.copy()
        D = list(set(A) | set(B) | set(C))

        E = Gene_list_pick (gene_list)[0]
        for i in range(36):
            coexp = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/Coexp_modules/adjacency_matrix_corr_module_M"+str(i+1)+".txt", sep=' ',low_memory=False)
            #select top 50% of correlations
            Final = coexp.sort_values(["corr"], ascending=False).reset_index(drop=True)
            Final2 = Final.head(int(0.5*len(Final)))
            Final2 = coexp.copy()    
            Final3 = Final2[Final2["row"].isin(D)&Final2["col"].isin(E)]
            Final4 = Final2[Final2["row"].isin(E)&Final2["col"].isin(D)]
            Final5 = Pd.concat([Final3,Final4])
            Final5["Status"] =Final5["row"].isin(E)
            Final5["ASD_gene"] = [d if e else c for d, c, e in zip(Final5.row,Final5.col,Final5.Status)]
            Final5["lncRNA"] = [c if e else d for d, c, e in zip(Final5.row,Final5.col,Final5.Status)]
    
            lncRNA = self.lncRNAKB_Final_updated.copy()
            Part_A = Final5[Final5["lncRNA"].isin(A)].reset_index(drop=True)
            Part_A = Part_A.rename(columns={"lncRNA":"lncRNAKB"})
            Part_A_new = (
                Part_A.merge(lncRNA, 
                             on=['lncRNAKB'],
                             how='left'))
            #Part_A_new = Part_A_new.rename(columns={"lncRNAKB":"lncRNA"})
    
            Part_B = Final5[Final5["lncRNA"].isin(B)].reset_index(drop=True)
            Part_B = Part_B.rename(columns={"lncRNA":"Lincpedia"})
            Part_B_new = (
                Part_B.merge(lncRNA,
                             on=['Lincpedia'],
                             how='left'))
            #Part_B_new = Part_B_new.rename(columns={"Lincpedia":"lncRNA"})
 
            Part_C = Final5[Final5["lncRNA"].isin(C)].reset_index(drop=True)
            Part_C = Part_C.rename(columns={"lncRNA":"GENE_symbol"})
            Part_C_new = (
                Part_C.merge(lncRNA,
                             on=['GENE_symbol'],
                             how='left'))
            #Part_C_new = Part_C_new.rename(columns={"GENE_symbol":"lncRNA"})

            Final6 = Pd.concat([Part_A_new,Part_B_new,Part_C_new]).reset_index(drop=True)
            Final7 = Final6.drop_duplicates().reset_index(drop=True)
            Final7.to_csv("/hpf/largeprojects/tcagstor/users/sghaffari/Coexp_modules/lncRNA_"+gene_list+"_coexp_M"+str(i+1)+".tsv", sep="\t")

        DFs=[]
        for i in range(36):
            DF = Pd.read_csv("/hpf/largeprojects/tcagstor/users/sghaffari/Coexp_modules/lncRNA_"+gene_list+"_coexp_M"+str(i+1)+".tsv", sep="\t")
            DFs.append(DF)
        master_df = Pd.concat(DFs)
        master_df = master_df.drop_duplicates().reset_index(drop=True)
        return master_df

        

