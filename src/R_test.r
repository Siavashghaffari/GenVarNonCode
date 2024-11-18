#!/usr/bin/Rscript --vanilla
df = read.table(file = 'lncRNA_Overlapping_CNVs_NO_EXPRESSION.tsv', sep = '\t', header = TRUE)
df[is.na(df)] <- 0

p_values <- function(x1,x2,n1,n2) {prop.test(c(x1, x2), c(n1, n2), 
                                               alternative="greater", conf.level = 0.8, correct=F)$p.value}

# P-values for total cases vs. control
df$pVal_siblings <- mapply(p_values, df$N_probands, df$N_unaffectedsibling, df$all_probands, df$all_siblings)
df$pVal_1000G <- mapply(p_values, df$N_probands, df$N_1000G, df$all_probands, df$all_siblings)

# P-values for male only cases vs. control
df$pVal_siblings_male <- mapply(p_values, df$N_probands_male, df$N_unaffected_male, df$male_probands, df$male_siblings)
df$pVal_1000G_male <- mapply(p_values, df$N_probands_male, df$N_Control_male, df$all_probands, df$all_siblings)

# q-values for total cases vs. control
df$FDR_siblings =  p.adjust(df$pVal_siblings, method = "BH")
df$FDR_1000G =  p.adjust(df$pVal_1000G, method = "BH")

# q-values for male only cases vs. control
df$FDR_siblings_male =  p.adjust(df$pVal_siblings_male, method = "BH")
df$FDR_1000G_male =  p.adjust(df$pVal_1000G_male, method = "BH")

write.csv(df,file="P-values_lncRNA_NO_Expression_CNVs.csv",row.names = FALSE)
