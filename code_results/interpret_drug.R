library(readr)
library(tidyverse)
library(ggplot2)



drugcell_all <- read_delim("data_newgenes_ge_ont_drfeat/drugcell_all.txt",
                           "\t", escape_double = FALSE, col_names = FALSE,
                           trim_ws = TRUE)
colnames(drugcell_all) = c("ccl", "smiles", "auc")

drug_transl_pubchem <- read_delim("drug_transl_pubchem.csv",
                                  "\t", escape_double = FALSE, trim_ws = TRUE)

drugcell_all = drugcell_all %>% left_join(drug_transl_pubchem, by = "smiles")
drugcell_all$idx = 1:nrow(drugcell_all)
drugcell_all$family <- sub(".+?_", "", drugcell_all$ccl)

if(F) {
    DeepLift <- read_csv("code_results/interpretation/DeepLift.csv")
    DeepLift$idx = 1:nrow(DeepLift)

    DeepLift_wo_MF <- read_csv("code_results/interpretation/DeepLift_movida_newdrfeat_all_wo_MF.csv")
    DeepLift_wo_MF$idx = 1:nrow(DeepLift_wo_MF)

    DeepLift_pcfp <- read_csv("code_results/interpretation/DeepLift_movida_newdrfeat_all_pcfp.csv")
    DeepLift_pcfp$idx = 1:nrow(DeepLift_pcfp)
}


df_interpret = DeepLift_pcfp



idx_sel = drugcell_all %>% subset(ccl == "201T_LUNG")
idx_sel = idx_sel$idx

int_sel = df_interpret[idx_sel,]




mean_abs_columns = apply(df_interpret[,-ncol(df_interpret)], 2, function(x) mean(abs(x)))
mean_abs_columns_sorted = sort(mean_abs_columns, decreasing = T)








#int_sel[, colSums(abs(int_sel) > 0.05) > 0]
colnames(int_sel)[colSums(abs(int_sel) > 0.15) > 0]
int_col_imp = colSums(abs(int_sel) > 0.15)
int_col_imp = int_col_imp[-length(int_col_imp)]
int_col_imp = int_col_imp[int_col_imp > 0]
nrow(int_sel)
sort(int_col_imp, decreasing = T)

#min(df_interpret[,-ncol(df_interpret)])
#max(df_interpret[,-ncol(df_interpret)])



data_long <- int_sel[, -ncol(int_sel)] %>%
    pivot_longer(colnames(int_sel)[-ncol(int_sel)]) %>%
    as.data.frame()
ggplot(data_long, aes(x = value)) +
    geom_density()








