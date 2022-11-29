##################################################################################################################
############## Imports
##################################################################################################################
library(readr)
library(tidyverse)
library(yaGST)
library(writexl)
library(readxl)

##################################################################################################################
############## Parameters
##################################################################################################################

model_name = "model1_falseweighted_80eps_300e_mmse"
max_val = 1.2

top_x = 50
auc_mid = 0.3
drug_rf = "doxorubicin"
family_ccl = "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"
read_csv_bool = F
scale_bool = F

drug_bool = drug_rf != ""
family_bool = family_ccl != ""

if(F) {
    rm(list=setdiff(ls(), "silent_locimpr_val"))
}



##################################################################################################################
############## Read of inputs
############## drug_transl_pubchem      traslation of drug name
############## gos_descr                traslation of go description
############## drugcell_all             all triplets ccl, cmp, AUC, with AUC_pred, family name, drug name
############## silent_val               silent score for each triplets-go
############## silent_sel               selent_val subsetted for a drug, if any
############## silent_sel_family        selent_sel subseted for a family, if any
##################################################################################################################

drug_transl_pubchem <- read_delim("drug_transl_pubchem.csv",
                                  "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(drug_transl_pubchem) <- c("cmp", "name")

gos_descr <- read_csv("code_results/go_descr.csv")
pred <- read_csv(paste0("launch/all_tests_newgenes/", model_name,"/Result2/drugcell.predict"),
                 col_names = FALSE)

drugcell_all <- read_delim("data_newgenes_ge_ont/drugcell_all.txt",
                           "\t", escape_double = FALSE, col_names = FALSE,
                           trim_ws = TRUE)
colnames(drugcell_all) = c("ccl", "cmp", "AUC")
drugcell_all$AUC_pred = pred$X1
drugcell_all$X = 1:nrow(drugcell_all)

drugcell_all$AUC[drugcell_all$AUC >= max_val] <- max_val - 0.000001
drugcell_all$AUC_pred[drugcell_all$AUC_pred < 0] <- 0
drugcell_all$AUC_pred[drugcell_all$AUC_pred >= max_val] <- max_val - 0.000001


drugcell_all$family <- sub(".+?_", "", drugcell_all$ccl)
drugcell_all = drugcell_all %>% left_join(drug_transl_pubchem, by = "cmp")

families = unique(drugcell_all$family)
drugs = unique(drugcell_all$name)

if (read_csv_bool) {
    file_names <- list.files(path="code_results/silent_locimpr_val", pattern=".csv")
    silent_locimpr_val <- do.call(rbind, lapply(paste0("code_results/silent_locimpr_val/", file_names), read.csv))
    colnames(silent_locimpr_val) = gsub("\\.", ":", colnames(silent_locimpr_val))
    silent_locimpr_val$X = 1:nrow(silent_locimpr_val)
}


if (drug_bool) {
    sample_sel <- drugcell_all %>% subset(name == drug_rf)
} else {
    sample_sel <- drugcell_all
}

silent_sel <- silent_locimpr_val %>% subset(X %in% sample_sel$X)



if (family_bool) {
    sample_sel_family = sample_sel %>% subset(family == family_ccl)
} else {
    sample_sel_family = sample_sel
}

silent_sel_family <- silent_sel %>% subset(X %in% sample_sel_family$X)




##################################################################################################################
############ most_important_go function
##############
############## rank_topGO_foreach_ccldrug             rank of top_x GO for each triplets in training
############## frenquency_topGO_foreach_ccldrug       frequency of top_x GO for each triplets in training
############## mean_topGO                             rank of top_x GO for mean silent score
##################################################################################################################
#81044
most_important_go <- function(silent_sel_family_sub, sc_bool, top_x, gos_descr, drug_rf, family_ccl, str_sub, str_sil = "silent_score") {

    tmp = as.data.frame(t(silent_sel_family_sub[,2:ncol(silent_sel_family_sub)]))
    tmp = tmp[row.names(tmp) != "GO:0008150", , drop=F]
    if (scale_bool) {
        tmp = scale(tmp)
    }

    ord = apply(tmp, 2, function(x) order(x, decreasing = T))
    if (!is.null(top_x)) {
        ord = ord[1:top_x, , drop=F]
    }
    rank_topGO_foreach_ccldrug = apply(ord, 2, function(x) rownames(tmp)[x])
    frenquency_topGO_foreach_ccldrug = sort(table(rank_topGO_foreach_ccldrug))



    sub_mean = as.data.frame(rowMeans(tmp))
    sub_mean$go = rownames(sub_mean)
    colnames(sub_mean) = c("silent_score", "go")
    rownames(sub_mean) = NULL

    res = sub_mean
    res <- res %>% inner_join(gos_descr, by = "go")
    res$silent_score <- as.double(res$silent_score)

    res2 <- res
    res2 <- within(res2, words <- gsub("(([A-Za-z1-9.,']+\\s){3})","\\1\n ", descr))

    if (is.null(top_x)) {
        p = res2 %>%
            arrange(desc(silent_score)) %>%
            ggplot(., aes(x = reorder(words, -silent_score) , y = silent_score)) +
            geom_bar(stat='identity', color = "cornflowerblue", fill = "cornflowerblue", alpha = 0.4) +
            theme_bw() +
            labs(title = paste('Mean', str_sil, 'for', str_sub ,'applying', drug_rf, "to", family_ccl), x = "GO", y = "Silent score") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15))
    } else {
        p = res2 %>%
            arrange(desc(silent_score)) %>%
            head(top_x) %>%
            ggplot(., aes(x = reorder(words, -silent_score) , y = silent_score)) +
            geom_bar(stat='identity', color = "cornflowerblue", fill = "cornflowerblue", alpha = 0.4) +
            theme_bw() +
            labs(title = paste('Mean', str_sil, 'for', str_sub ,'applying', drug_rf, "to", family_ccl), x = "GO", y = "Silent score") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15))
    }
    #print(p)

    res = res[order(-res2$silent_score),]

    if (is.null(top_x)) {
        mean_topGO = res
    } else {
        mean_topGO = res[1:top_x,]
    }
    list("mean_topGO" = mean_topGO, "freq_go" = frenquency_topGO_foreach_ccldrug,
         "rank_go" = rank_topGO_foreach_ccldrug, "plot" = p)
}


##################################################################################################################
############## Most important GO for all considered drug-family
##################################################################################################################
sample_sel2 = sample_sel_family#[1,]
silent_locimpr_sel = silent_locimpr_val %>% subset(X %in% sample_sel2$X)


#silent_locimpr_sel = silent_locimpr_val[129, ]


res_tmp = most_important_go(silent_locimpr_sel, scale_bool, 200, gos_descr, drug_rf, family_ccl, "all", "silent_locimp")
mean_topGO_all = res_tmp$mean_topGO
frenquency_topGO_foreach_ccldrug_all = res_tmp$freq_go
rank_topGO_foreach_ccldrug_all = res_tmp$rank_go

p = res_tmp$plot
aspect_ratio = 3
#print(p)
ggsave("./code_results/silent_locimpr_score.pdf", p, height = 7 , width = 10 * aspect_ratio)




##################################################################################################################
############## mwwGST
##################################################################################################################


comp_gsea = function(silent_locimpr_compl, drugcell_all, gos_descr, drug_cons, auc_mid, family_ccl = "", thr_filter = 20, head_x = 50) {
    #print(drug_cons)
    silent_locimpr_compl = cbind(silent_locimpr_compl, drugcell_all[, c(1,7,3,4,6)])

    silent_locimpr_sel = silent_locimpr_compl %>% subset(name == drug_cons)

    if (family_ccl != "") {
        silent_locimpr_sel = silent_locimpr_sel %>% subset(family == family_ccl)
    }

    if (nrow(silent_locimpr_sel) < 50) {
        NULL
    } else {
        silent_locimpr_sel_ord = silent_locimpr_sel %>% dplyr::select("ccl", "name", "AUC", "AUC_pred", "family", "X", everything())
        silent_locimpr_sel_tr = as.data.frame(t(silent_locimpr_sel_ord))
        silent_locimpr_sel_tr = silent_locimpr_sel_tr[row.names(silent_locimpr_sel_tr) != "GO:0008150", , drop=F]

        silent_locimpr_sel_go = silent_locimpr_sel_tr[7:nrow(silent_locimpr_sel_tr), , drop=FALSE]
        silent_locimpr_sel_go = mutate_all(silent_locimpr_sel_go, function(x) as.numeric(as.character(x)))
        silent_locimpr_sel_go$go = rownames(silent_locimpr_sel_go)
        rownames(silent_locimpr_sel_go) = NULL
        silent_locimpr_sel_go = silent_locimpr_sel_go %>% left_join(gos_descr, by = "go")
        silent_locimpr_sel_go = silent_locimpr_sel_go %>% dplyr::select("go", "descr", everything())

        silent_locimpr_sel_descr = silent_locimpr_sel_tr[1:6, , drop=FALSE]


        ord = order(silent_locimpr_sel_descr["AUC_pred",])
        ccl_ord = unlist(silent_locimpr_sel_descr["ccl", ord])

        silent_locimpr_sel_go = silent_locimpr_sel_go[, c("go", "descr", colnames(silent_locimpr_sel_descr)[ord])]

        silent_locimpr_sel_descr = silent_locimpr_sel_descr[, ord, drop=FALSE]

        #View(silent_locimpr_sel_go[, 1:12])
        #View(silent_locimpr_sel_descr[, 1:10])

        silent_locimpr_sel_descr_tr = as.data.frame(t(silent_locimpr_sel_descr))
        silent_locimpr_sel_descr_tr$AUC = as.double(silent_locimpr_sel_descr_tr$AUC)
        silent_locimpr_sel_descr_tr$AUC_pred = as.double(silent_locimpr_sel_descr_tr$AUC_pred)
        sn = sum(silent_locimpr_sel_descr_tr$AUC_pred < auc_mid)
        #sn = floor(nrow(silent_locimpr_sel_descr_tr) / 10)
        #print(nrow(silent_locimpr_sel_descr_tr))
        #print(sn)
        #if ((sn >= max(floor(nrow(silent_locimpr_sel_descr_tr)/10), 1) | sn > 20) & nrow(silent_locimpr_sel_descr_tr) > sn + 1) {
        if (sn >= 5 & nrow(silent_locimpr_sel_descr_tr) > sn + 1) {
            res_gsea = apply(silent_locimpr_sel_go[, 3:ncol(silent_locimpr_sel_go)], 1, mwwGST, colnames(silent_locimpr_sel_go)[3:(sn+2)], 1)
            nes_pval = res_gsea %>% map_dfr(function(x) x[c("nes", "p.value")])
            nes_pval$go = silent_locimpr_sel_go$go
            nes_pval$descr = silent_locimpr_sel_go$descr

            nes_pval$p.value.adj = p.adjust(nes_pval$p.value, method = "fdr")

            nes_pval_sel = nes_pval %>% filter(p.value.adj < 0.05) %>% arrange(desc(nes))
            if (!(is.null(head_x))) {
                nes_pval_sel = nes_pval_sel %>% head(head_x)
            }
            nes_pval_sel$drug = drug_cons

            nes_pval_sel
        }
    }


}


silent_locimpr_completo = silent_locimpr_val

# HAEMATOPOIETIC_AND_LYMPHOID_TISSUE
nes_pval = comp_gsea(silent_locimpr_completo, drugcell_all, gos_descr, "doxorubicin", auc_mid, "LUNG", 15)


##################################################################################################################
############## Fisher read inputs
##################################################################################################################

cell2mutation <- read_csv("data_newgenes_ge_ont/cell2mutation.txt",
                          col_names = FALSE)
cell2amp <- read_csv("data_newgenes_ge_ont/cell2amp.txt",
                          col_names = FALSE)
cell2del <- read_csv("data_newgenes_ge_ont/cell2del.txt",
                          col_names = FALSE)
cell2ind <- read_delim("data_newgenes_ge_ont/cell2ind.txt",
                       "\t", escape_double = FALSE, col_names = FALSE,
                       trim_ws = TRUE)
gene2ind2mut <- read_delim("data_newgenes_ge_ont/gene2ind2mut.txt",
                           "\t", escape_double = FALSE, col_names = FALSE,
                           trim_ws = TRUE)
gene2ind2amp <- read_delim("data_newgenes_ge_ont/gene2ind2amp.txt",
                           "\t", escape_double = FALSE, col_names = FALSE,
                           trim_ws = TRUE)
gene2ind2del <- read_delim("data_newgenes_ge_ont/gene2ind2del.txt",
                           "\t", escape_double = FALSE, col_names = FALSE,
                           trim_ws = TRUE)
ont_mut <- read_delim("data_newgenes_ge_ont/ont_mut.txt",
                      "\t", escape_double = FALSE, col_names = FALSE,
                      trim_ws = TRUE)
ont_amp <- read_delim("data_newgenes_ge_ont/ont_amp.txt",
                      "\t", escape_double = FALSE, col_names = FALSE,
                      trim_ws = TRUE)
ont_del <- read_delim("data_newgenes_ge_ont/ont_del.txt",
                      "\t", escape_double = FALSE, col_names = FALSE,
                      trim_ws = TRUE)

ont_amp = unique(ont_amp)
ont_del = unique(ont_del)

cell2mut = as.data.frame(cell2mutation)
colnames(cell2mut) = gene2ind2mut$X2
rownames(cell2mut) = cell2ind$X2

cell2amp = as.data.frame(cell2amp)
colnames(cell2amp) = gene2ind2amp$X2
rownames(cell2amp) = cell2ind$X2

cell2del = as.data.frame(cell2del)
colnames(cell2del) = gene2ind2del$X2
rownames(cell2del) = cell2ind$X2

colnames(cell2ind) = c("X1", "ccl")
cell2ind$family <- sub(".+?_", "", cell2ind$ccl)

##################################################################################################################
############## Fisher
##################################################################################################################


fisher_gene = function(drug_sel, fam_sel, gene_sel, cell2ind, cell2x, drugcell_all) {
    ccl_sel_fam = cell2ind
    if (fam_sel != "") {
        ccl_sel_fam = cell2ind %>% subset(family == fam_sel)
        drugcell_sel = drugcell_all %>% subset(family == fam_sel)
    }
    x_sel_fam = cell2x[ccl_sel_fam$ccl, gene_sel]
    ccl_sel_fam$gene_x = x_sel_fam

    #drugcell_sel = drugcell_all %>% subset(family == fam_sel & name == drug_sel)
    drugcell_sel = drugcell_all %>% subset(name == drug_sel)
    drugcell_sel = drugcell_sel[, c("ccl", "AUC_pred")]
    drugcell_sel = unique(drugcell_sel)
    drugcell_sel$thr = ifelse(drugcell_sel$AUC_pred < 0.3, 1, 0)

    drugcell_sel = drugcell_sel %>% left_join(ccl_sel_fam, by = "ccl")

    #drugcell_sel[,c("thr", "gene_x")]
    #count_thr_x = count(drugcell_sel, thr, gene_x, .drop = F)
    count_thr_x = table(drugcell_sel$thr, drugcell_sel$gene_x)
    if (nrow(count_thr_x) < 2 | ncol(count_thr_x) < 2) {
        NULL
    } else {
        fisher_res = fisher.test(count_thr_x)
        list("count_x" = count_thr_x, "pval" = fisher_res$p.value)
    }
}

fisher_gene_filtered = function(drug_sel, fam_sel, gene_sel, cell2ind, cell2x, drugcell_all, count_min = 10) {
    fisher_res_list = fisher_gene(drug_sel, fam_sel, gene_sel, cell2ind, cell2x, drugcell_all)
    if (is.null(fisher_res_list)) {
        NULL
    } else {
        pval_fisher = fisher_res_list$pval
        count_x = fisher_res_list$count_x
        n_x = sum(count_x[,2])
        if (pval_fisher < 0.05) { #& n_x >= floor(sum(count_x)/10)) {
            #fisher_res = rbind(fisher_res, c(go, gene_sel, pval_fisher))
            #count_mut_res = rbind(count_mut_res, count_mut)
            gene_sel
        } else {
            NULL
        }
    }

}

get_important_df = function(imp_df) {

    df_res = data.frame()
    for (df in imp_df) {
        if (! is.null(df)) {
            df_res = rbind(df_res, df)
        }
    }

    rownames(df_res) = 1:nrow(df_res)
    colnames(df_res) = c("drug", "family", "go", "descr", "nes", "genes_mut", "genes_amp", "genes_del")
    df_res
}

silent_locimpr_completo = silent_locimpr_val




important_go_gene_in_drugfamily = c()
cores = 80
doMC::registerDoMC(cores)
# families
important_all <- foreach (i=1:length(drugs)) %dopar% {
#}


    important_go_gene_in_drugfamily_tmp = c()
#for( i in length(drugs)) {

    drug_sel = drugs[i]
    print(paste(i, "/", length(drugs), drug_sel))


    for (family_sel in families) {
        print(paste(drug_sel, family_sel))
        nes_pval = comp_gsea(silent_locimpr_completo, drugcell_all, gos_descr, drug_sel, auc_mid, family_sel, 15)
        if (is.null(nes_pval)) {
            next
        }
        nes_pval_all_unlist = bind_rows(nes_pval, .id = "column_label")
        nes_pval_all_unlist = nes_pval_all_unlist[, -1]
        nes_pval_all_unlist = nes_pval_all_unlist[, c("go", "descr", "nes")]

        for (go in nes_pval_all_unlist$go) {

            gene_in_go = (ont_mut %>% filter(X1 == go))$X2
            genes_mut = c()
            if (length(gene_in_go) > 0) {
                for (gene_sel in gene_in_go) {
                    genes_filtered = fisher_gene_filtered(drug_sel, family_sel, gene_sel, cell2ind, cell2mut,
                                                          drugcell_all, count_min = 10)
                    if ( !is.null(genes_filtered) ) {
                        genes_mut = c(genes_mut, genes_filtered)
                    }
                }
            }


            gene_in_go = (ont_amp %>% filter(X1 == go))$X2
            genes_amp = c()
            if (length(gene_in_go) > 0) {
                for (gene_sel in gene_in_go) {
                    genes_filtered = fisher_gene_filtered(drug_sel, family_sel, gene_sel, cell2ind, cell2amp,
                                                          drugcell_all, count_min = 10)
                    if ( !is.null(genes_filtered) ) {
                        genes_amp = c(genes_amp, genes_filtered)
                    }
                }
            }


            gene_in_go = (ont_del %>% filter(X1 == go))$X2
            genes_del = c()
            if (length(gene_in_go) > 0) {
                for (gene_sel in gene_in_go) {
                    genes_filtered = fisher_gene_filtered(drug_sel, family_sel, gene_sel, cell2ind, cell2del,
                                                          drugcell_all, count_min = 10)
                    if ( !is.null(genes_filtered) ) {
                        genes_del = c(genes_del, genes_filtered)
                    }

                }
            }

            nes_pval_onego = nes_pval_all_unlist[nes_pval_all_unlist$go == go,]

            genes_mut = paste(genes_mut, collapse = ",")
            genes_amp = paste(genes_amp, collapse = ",")
            genes_del = paste(genes_del, collapse = ",")

            go_genes_mutampdel = c(drug_sel, family_sel, go, nes_pval_onego$descr, nes_pval_onego$nes,
                                   genes_mut, genes_amp, genes_del)
            important_go_gene_in_drugfamily_tmp = rbind(important_go_gene_in_drugfamily_tmp, go_genes_mutampdel)
        }


    }

    important_go_gene_in_drugfamily_tmp
}
# all
important_all <- foreach (i=1:length(drugs)) %dopar% {
    #}


    important_go_gene_in_drugfamily_tmp = c()
    #for( i in length(drugs)) {

    drug_sel = drugs[i]
    print(paste(i, "/", length(drugs), drug_sel))


    #for (family_sel in families) {
    family_sel = ""
    print(paste(drug_sel))#, family_sel))
    nes_pval = comp_gsea(silent_locimpr_completo, drugcell_all, gos_descr, drug_sel, auc_mid, family_sel, 15)
    if (!is.null(nes_pval)) {

        nes_pval_all_unlist = bind_rows(nes_pval, .id = "column_label")
        nes_pval_all_unlist = nes_pval_all_unlist[, -1]
        nes_pval_all_unlist = nes_pval_all_unlist[, c("go", "descr", "nes")]

        for (go in nes_pval_all_unlist$go) {

            gene_in_go = (ont_mut %>% filter(X1 == go))$X2
            genes_mut = c()
            if (length(gene_in_go) > 0) {
                for (gene_sel in gene_in_go) {
                    genes_filtered = fisher_gene_filtered(drug_sel, family_sel, gene_sel, cell2ind, cell2mut,
                                                          drugcell_all, count_min = 10)
                    if ( !is.null(genes_filtered) ) {
                        genes_mut = c(genes_mut, genes_filtered)
                    }
                }
            }


            gene_in_go = (ont_amp %>% filter(X1 == go))$X2
            genes_amp = c()
            if (length(gene_in_go) > 0) {
                for (gene_sel in gene_in_go) {
                    genes_filtered = fisher_gene_filtered(drug_sel, family_sel, gene_sel, cell2ind, cell2amp,
                                                          drugcell_all, count_min = 10)
                    if ( !is.null(genes_filtered) ) {
                        genes_amp = c(genes_amp, genes_filtered)
                    }
                }
            }


            gene_in_go = (ont_del %>% filter(X1 == go))$X2
            genes_del = c()
            if (length(gene_in_go) > 0) {
                for (gene_sel in gene_in_go) {
                    genes_filtered = fisher_gene_filtered(drug_sel, family_sel, gene_sel, cell2ind, cell2del,
                                                          drugcell_all, count_min = 10)
                    if ( !is.null(genes_filtered) ) {
                        genes_del = c(genes_del, genes_filtered)
                    }

                }
            }

            nes_pval_onego = nes_pval_all_unlist[nes_pval_all_unlist$go == go,]

            genes_mut = paste(genes_mut, collapse = ",")
            genes_amp = paste(genes_amp, collapse = ",")
            genes_del = paste(genes_del, collapse = ",")

            family_sel = "all"
            go_genes_mutampdel = c(drug_sel, family_sel, go, nes_pval_onego$descr, nes_pval_onego$nes,
                                   genes_mut, genes_amp, genes_del)
            important_go_gene_in_drugfamily_tmp = rbind(important_go_gene_in_drugfamily_tmp, go_genes_mutampdel)
        }
            #}

    }

    important_go_gene_in_drugfamily_tmp
}


important_go_gene_in_drugfamily = get_important_df(important_all)
View(important_go_gene_in_drugfamily)
write_csv(important_go_gene_in_drugfamily, "./code_results/important_go_gene_in_drugfamily.csv")



##################################################################################################################
############## Fisher heatmap
##################################################################################################################


important_sel = important_go_gene_in_drugfamily[, c("drug", "family", "descr")]
important_sel = important_sel %>% subset(family == "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE")
important_sel = unique(important_sel)
important_sel = important_sel[, c("drug", "descr")]
important_sel["val"] = 1
important_sel_wid = as.data.frame(important_sel %>% pivot_wider(names_from = descr, values_from = val, values_fill = 0))
rownames(important_sel_wid) = substr(important_sel_wid$drug, 0, 20)
important_sel_wid = important_sel_wid[, -1]
important_sel_wid = important_sel_wid[, colSums(important_sel_wid) > 3]
colnames(important_sel_wid) = substr(colnames(important_sel_wid), 0, 20)

d = dist(important_sel_wid, method="binary")
fit = hclust(d, method="ward.D2")
plot(fit, cex=1)
heatmaply(as.matrix(d))
heatmaply(as.matrix(important_sel_wid), showticklabels = T)












cores = 80
chunks = split(drugs, ceiling(seq_along(drugs)/20))
doMC::registerDoMC(cores)
nes_pval_all <- foreach (i=1:length(chunks)) %dopar% {
    chunk = chunks[i]
    chunk = chunk[[1]]
    lapply(chunk, function(x) comp_gsea(silent_locimpr_completo, drugcell_all, gos_descr,
                                        x, auc_mid, "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
                                        head_x = NULL))
}

nes_pval_all_unlist = bind_rows(nes_pval_all, .id = "column_label")
nes_pval_all_unlist = nes_pval_all_unlist[, -1]
#write_csv(nes_pval_all_unlist, "./code_results/nes_pval.csv")

nes_pval_sel = nes_pval_all_unlist[, c("drug", "descr", "nes")]
nes_pval_sel = unique(nes_pval_sel)
nes_pval_sel_wid = as.data.frame(nes_pval_sel %>%
                                      pivot_wider(names_from = descr, values_from = nes, values_fill = 0))

rownames(nes_pval_sel_wid) = substr(nes_pval_sel_wid$drug, 0, 20)
nes_pval_sel_wid = nes_pval_sel_wid[, -1]
nes_pval_sel_wid = nes_pval_sel_wid[, colSums(nes_pval_sel_wid) > 3]
colnames(nes_pval_sel_wid) = substr(colnames(nes_pval_sel_wid), 0, 20)

d = dist(nes_pval_sel_wid, method="euclidian")
fit = hclust(d, method="ward.D2")
plot(fit, cex=1)
heatmaply(as.matrix(d))
heatmaply(as.matrix(nes_pval_sel_wid), showticklabels = T)

##################################################################################################################
########### Drug Target and Drugs Synergy DB
##################################################################################################################


drug_target <- read_csv("code_results/ttd/drug_target.csv")
drug_target = unique(drug_target)

breast_anchor_combo <- read_csv("code_results/ttd/breast_anchor_combo.csv")
colon_anchor_combo <- read_csv("code_results/ttd/colon_anchor_combo.csv")
pancreas_anchor_combo <- read_csv("code_results/ttd/pancreas_anchor_combo.csv")




ont_mut <- read_delim("data_newgenes_ge_ont_drfeat/ont_mut.txt",
                      "\t", escape_double = FALSE, col_names = FALSE,
                      trim_ws = TRUE)
ont_del <- read_delim("data_newgenes_ge_ont_drfeat/ont_del.txt",
                      "\t", escape_double = FALSE, col_names = FALSE,
                      trim_ws = TRUE)
ont_amp <- read_delim("data_newgenes_ge_ont_drfeat/ont_amp.txt",
                      "\t", escape_double = FALSE, col_names = FALSE,
                      trim_ws = TRUE)

colnames(ont_mut) = c("GO", "gene")
colnames(ont_del) = c("GO", "gene")
colnames(ont_amp) = c("GO", "gene")

ont_gene = rbind(ont_mut, ont_del, ont_amp)
ont_gene = unique(ont_gene)


synergy_drugs = rbind(breast_anchor_combo, colon_anchor_combo, pancreas_anchor_combo)
synergy_drugs$`Cell Line name` = gsub("-", "", synergy_drugs$`Cell Line name`)
synergy_drugs$`Cell Line name` = toupper(synergy_drugs$`Cell Line name`)

cell_lines = unique(drugcell_all$ccl)
cell_lines = gsub("\\_.*","", cell_lines)



screened_drugs_syn = synergy_drugs[, c("Anchor Name", "Anchor Target")]
colnames(screened_drugs_syn) = c("drug", "target")
screened_drugs_syn_tmp = synergy_drugs[, c("Library Name", "library Target")]
colnames(screened_drugs_syn_tmp) = c("drug", "target")
screened_drugs_syn = rbind(screened_drugs_syn, screened_drugs_syn_tmp)
screened_drugs_syn = unique(screened_drugs_syn)
rm(screened_drugs_syn_tmp)




drug_target = drug_target[drug_target$drug %in% screened_drugs_syn$drug, ]
drug_target = separate_rows(drug_target, "gene", sep = "-")


screened_drugs_syn[screened_drugs_syn$drug == "Vorinostat", "target"] = "HDAC"
screened_drugs_syn = separate_rows(screened_drugs_syn, "target", sep = "\\|")
screened_drugs_syn = separate_rows(screened_drugs_syn, "target", sep = ",")
screened_drugs_syn$target = str_trim(screened_drugs_syn$target)

colnames(screened_drugs_syn)[2] = "gene"


drug_target_all = rbind(drug_target[, 1:2], screened_drugs_syn)
drug_target_all = unique(drug_target_all)





synergy_drugs_feat = unique(synergy_drugs[, c("Cell Line name", "Anchor Name",
                                              "Library Name", "Synergy?", "Delta Xmid", "Delta Emax")])

synergy_drugs_feat = synergy_drugs_feat %>%
    rowwise() %>%
    mutate(drug1 = sort(c(`Anchor Name`, `Library Name`))[1],
           drug2 = sort(c(`Anchor Name`, `Library Name`))[2])

synergy_drugs_feat = synergy_drugs_feat[, c(1, 7, 8, 4, 5, 6)]
colnames(synergy_drugs_feat)[1] = "ccl"

synergy_drugs_mng = synergy_drugs_feat %>% group_by(ccl, drug1, drug2) %>%
    summarise(s = min(sum(`Synergy?`), 1),
              max_D_Xmid = max(`Delta Xmid`),
              max_D_Emax = max(`Delta Emax`))


##################################################################################################################
########### Predicting synergy by ccl
##################################################################################################################






tmp = synergy_drugs_mng
tmp1 = tmp[, 1:2]
tmp2 = tmp[, c(1,3)]
colnames(tmp1)[2] = "drug"
colnames(tmp2)[2] = "drug"
ass_ccl_dr = rbind(tmp1, tmp2)
ass_ccl_dr = unique(ass_ccl_dr)

rm(tmp, tmp1, tmp2)

ass_ccl_dr_sel = ass_ccl_dr %>% subset(ccl %in% cell_lines)
ass_ccl_dr_sel = ass_ccl_dr_sel %>% subset(tolower(drug) %in% tolower(drugs))



doMC::registerDoMC(80)


res_par = foreach::foreach(i = 1:nrow(ass_ccl_dr_sel)) %dopar% {
    ccl_cons = ass_ccl_dr_sel[[i, "ccl"]]
    drug_cons = ass_ccl_dr_sel[[i, "drug"]]






    sample_sel2 = drugcell_all
    sample_sel2 = sample_sel2 %>% filter(startsWith(ccl, ccl_cons), tolower(name) == tolower(drug_cons)) %>% head(1)

    if (nrow(sample_sel2) == 0) {
        c(drug_cons, ccl_cons, NA, NA)

    } else {
        silent_locimpr_sel = silent_locimpr_val %>% subset(X %in% sample_sel2$X)



        res_tmp = most_important_go(silent_locimpr_sel, scale_bool, NULL, gos_descr, drug_cons, family_ccl, "all", "silent_locimp")
        mean_topGO_all = res_tmp$mean_topGO


        mean_topGO_all_high = mean_topGO_all %>% head(5)

        genes_sel = unique((ont_gene %>% filter(GO %in% mean_topGO_all_high$go))$gene)


        drug_target_sel_high = drug_target_all %>% subset(gene %in% genes_sel | gene %in% gsub("\\d+$", "", genes_sel))

        drug_target_sel_high$drug_c = drug_cons
        drug_target_sel_high$ccl = ccl_cons
        drug_target_sel_high$s = 1



        mean_topGO_all_low = mean_topGO_all %>% tail(5)

        genes_sel = unique((ont_gene %>% filter(GO %in% mean_topGO_all_low$go))$gene)


        drug_target_sel_low = drug_target_all %>% subset(gene %in% genes_sel | gene %in% gsub("\\d+$", "", genes_sel))

        drug_target_sel_low$drug_c = drug_cons
        drug_target_sel_low$ccl = ccl_cons
        drug_target_sel_low$s = 0
        #drug_target_sel = drug_target_sel[ , c("ccl", "drug_cons", "drug", "gene")]

        drug_target_sel = rbind(drug_target_sel_high, drug_target_sel_low)


        drug_target_sel = drug_target_sel %>%
            rowwise() %>%
            mutate(drug1 = sort(c(drug, drug_c))[1],
                   drug2 = sort(c(drug, drug_c))[2])

        drug_target_sel = drug_target_sel[, c("ccl", "drug1", "drug2", "s", "gene")]


        synergy_sel = synergy_drugs_mng
        synergy_sel = synergy_sel %>% filter(ccl == ccl_cons)
        synergy_sel = synergy_sel %>% filter(drug1 == drug_cons | drug2 == drug_cons)
        synergy_sel_wo_na = drug_target_sel %>% full_join(synergy_sel, by = c("ccl", "drug1", "drug2"))


        synergy_sel_wo_na = synergy_sel_wo_na[!is.na(synergy_sel_wo_na$s.x),]
        synergy_sel_wo_na = synergy_sel_wo_na[!is.na(synergy_sel_wo_na$s.y),]
        colnames(synergy_sel_wo_na)[4] = "syn_pred"
        colnames(synergy_sel_wo_na)[6] = "syn_truth"

        synergy_sel_wo_na = synergy_sel_wo_na[,c(1,2,3,4,6)]
        synergy_sel_wo_na = unique(synergy_sel_wo_na)
        t = table(synergy_sel_wo_na[, c("syn_pred", "syn_truth")])

        if (nrow(t) != 2 | ncol(t) != 2) {
            c(drug_cons, ccl_cons, NA, NA)

        } else {
            n_1 = sum(synergy_sel$s == 1)
            n_0 = sum(synergy_sel$s == 0)

            p = n_1/(n_1 + n_0)

            t2 = rbind(c(n_0, n_1), t[2,])

            bin_pval = binom.test(t[2,2], sum(t[2,]), p)$p.value
            fish_pval = fisher.test(t2)$p.value

            enr = (t2[2,2]/sum(t2[2,])) / (t2[1,2]/sum(t2[1,]))

            c(drug_cons, ccl_cons, fish_pval, bin_pval, enr)
        }
    }

}


synergy_res_ccl_pval = do.call(rbind.data.frame, res_par)

colnames(synergy_res_ccl_pval) = c("drug", "ccl", "fisher_pval", "binomial_pval", "enr")
synergy_res_ccl_pval$fisher_pval = as.numeric(synergy_res_ccl_pval$fisher_pval)
synergy_res_ccl_pval$binomial_pval = as.numeric(synergy_res_ccl_pval$binomial_pval)
synergy_res_ccl_pval$enr = as.numeric(synergy_res_ccl_pval$enr)
#synergy_res_ccl_pval


write_csv(synergy_res_ccl_pval, "code_results/synergy/synergy_res_ccl_pval.csv")


synergy_res_ccl_pval_filtered = synergy_res_ccl_pval %>% filter(fisher_pval < 0.05 | binomial_pval < 0.05)








synergy_res_ccl_pval_wo_na = synergy_res_ccl_pval %>% filter(!is.na(fisher_pval) & !is.na(binomial_pval))

synergy_res_ccl_pval_wo_na$prob = 0.0

for (i in 1:nrow(synergy_res_ccl_pval_wo_na)) {
    synergy_sel = synergy_drugs_mng
    synergy_sel = synergy_sel %>% filter(ccl == synergy_res_ccl_pval_wo_na[i,"ccl"])
    synergy_sel = synergy_sel %>% filter(drug1 == synergy_res_ccl_pval_wo_na[i,"drug"] |
                                             drug2 == synergy_res_ccl_pval_wo_na[i,"drug"])
    n_1 = sum(synergy_sel$s == 1)
    n_0 = sum(synergy_sel$s == 0)

    p = n_1/(n_1 + n_0)
    synergy_res_ccl_pval_wo_na[i,"prob"] = p
}


dfpl = synergy_res_ccl_pval_wo_na
dfpl = dfpl %>% filter(!enr == 0, !fisher_pval == 1)

ggplot(data=dfpl, aes(x = log10(enr), y = -log10(fisher_pval))) +
    geom_point() +
    theme_minimal() +
    geom_hline(yintercept = -log10(0.05))

ggplot(data = dfpl, aes(x = log10(enr), y = -log10(fisher_pval))) +
    geom_point(alpha = 0.3) +
    geom_point(data = (dfpl %>% filter(-log10(fisher_pval) > -log10(0.05))),
               aes(x = log10(enr), y = -log10(fisher_pval)),
               color = 'green',
               size = 3) +
    geom_text(data = (dfpl %>% filter(-log10(fisher_pval) > -log10(0.05))),
              aes(x = log10(enr), y = -log10(fisher_pval), label = paste0(drug, " ", ccl)),
              position=position_jitter(width=0.1,height=0.1)) +
    theme_minimal()

ggsave("./code_results/synergy/volcano_plot.pdf", device = "pdf")

##################################################################################################################
########### Predicting synergy by tumor
##################################################################################################################


tmp = synergy_drugs_mng
tmp1 = tmp[, 1:2]
tmp2 = tmp[, c(1,3)]
colnames(tmp1)[2] = "drug"
colnames(tmp2)[2] = "drug"
ass_ccl_dr = rbind(tmp1, tmp2)
ass_ccl_dr = unique(ass_ccl_dr)

rm(tmp, tmp1, tmp2)



tmp = synergy_drugs[, c(1,3)]
colnames(tmp) = c("ccl", "family")
tmp = unique(tmp)
tmp$family = toupper(tmp$family)
tmp$family = gsub(" ", "_", tmp$family)

ass_ccl_dr_fam = ass_ccl_dr %>% left_join(tmp, by = "ccl")


ass_ccl_dr_fam = unique(ass_ccl_dr_fam[, c(2,3)])


ass_ccl_dr_fam_sel = ass_ccl_dr_fam %>% subset(tolower(drug) %in% tolower(drugs))


synergy_drugs_fam = synergy_drugs_mng %>% left_join(tmp, by = "ccl")
synergy_drugs_fam = unique(synergy_drugs_fam[, c(7,2,3,4)])


synergy_drugs_fam = synergy_drugs_fam %>% group_by(family, drug1, drug2) %>%
    summarise(s = min(sum(s), 1))





doMC::registerDoMC(6)


res_par = foreach::foreach(i = 1:nrow(ass_ccl_dr_fam_sel)) %dopar% {

    fam_cons = ass_ccl_dr_fam_sel[[i, "family"]]
    drug_cons = ass_ccl_dr_fam_sel[[i, "drug"]]

    drugcell_sel = drugcell_all %>% filter( family == fam_cons , tolower(name) == tolower(drug_cons))
    auc_mid = (floor(min(drugcell_sel$AUC_pred) * 10) + 3) / 10



    silent_locimpr_completo = silent_locimpr_val
    nes_pval = comp_gsea(silent_locimpr_completo, drugcell_all, gos_descr, drug_cons,
                         auc_mid, fam_cons, 5, head_x = NULL, to_low = T, th_num_sampl = F,
                         filtering = F)



    nes_pval_high = nes_pval %>% head(5)


    mean_topGO_all = nes_pval_high

    genes_sel = unique((ont_gene %>% filter(GO %in% mean_topGO_all$go))$gene)


    drug_target_sel_high = drug_target_all %>% subset(gene %in% genes_sel | gene %in% gsub("\\d+$", "", genes_sel))

    drug_target_sel_high$drug_c = drug_cons
    drug_target_sel_high$family = fam_cons
    drug_target_sel_high$s = 1



    nes_pval_low = nes_pval %>% tail(5)


    mean_topGO_all = nes_pval_low

    genes_sel = unique((ont_gene %>% filter(GO %in% mean_topGO_all$go))$gene)


    drug_target_sel_low = drug_target_all %>% subset(gene %in% genes_sel | gene %in% gsub("\\d+$", "", genes_sel))

    drug_target_sel_low$drug_c = drug_cons
    drug_target_sel_low$family = fam_cons
    drug_target_sel_low$s = 0
    #drug_target_sel = drug_target_sel[ , c("ccl", "drug_cons", "drug", "gene")]

    drug_target_sel = rbind(drug_target_sel_high, drug_target_sel_low)


    drug_target_sel = drug_target_sel %>%
        rowwise() %>%
        mutate(drug1 = sort(c(drug, drug_c))[1],
               drug2 = sort(c(drug, drug_c))[2])

    drug_target_sel = drug_target_sel[, c("family", "drug1", "drug2", "s", "gene")]


    synergy_sel = synergy_drugs_fam
    synergy_sel = synergy_sel %>% filter(family == fam_cons)
    synergy_sel = synergy_sel %>% filter(drug1 == drug_cons | drug2 == drug_cons)
    synergy_sel_wo_na = drug_target_sel %>% full_join(synergy_sel, by = c("family", "drug1", "drug2"))


    synergy_sel_wo_na = synergy_sel_wo_na[!is.na(synergy_sel_wo_na$s.x),]
    synergy_sel_wo_na = synergy_sel_wo_na[!is.na(synergy_sel_wo_na$s.y),]
    colnames(synergy_sel_wo_na)[4] = "syn_pred"
    colnames(synergy_sel_wo_na)[6] = "syn_truth"

    synergy_sel_wo_na = synergy_sel_wo_na[,c(1,2,3,4,6)]
    synergy_sel_wo_na = unique(synergy_sel_wo_na)
    t = table(synergy_sel_wo_na[, c("syn_pred", "syn_truth")])

    if (nrow(t) != 2 | ncol(t) != 2) {
        c(drug_cons, fam_cons, NA, NA)

    } else {
        n_1 = sum(synergy_sel$s == 1)
        n_0 = sum(synergy_sel$s == 0)

        p = n_1/(n_1 + n_0)

        t2 = rbind(c(n_0, n_1), t[2,])

        bin_pval = binom.test(t[2,2], sum(t[2,]), p)$p.value
        fish_pval = fisher.test(t2)$p.value

        enr = (t2[2,2]/sum(t2[2,])) / (t2[1,2]/sum(t2[1,]))

        c(drug_cons, fam_cons, fish_pval, bin_pval, enr)
    }

}

synergy_res_fam_pval = do.call(rbind.data.frame, res_par)

colnames(synergy_res_fam_pval) = c("drug", "family", "fisher_pval", "binomial_pval", "enr")

synergy_res_fam_pval


write_csv(synergy_res_fam_pval, "code_results/synergy/synergy_res_fam_pval.csv")

synergy_res_fam_pval_filtered = synergy_res_fam_pval %>% filter(fisher_pval < 0.05 | binomial_pval < 0.05)








rrr = rbind(breast_anchor_combo, colon_anchor_combo, pancreas_anchor_combo)

rrr = rrr %>%
    rowwise() %>%
    mutate(drug1 = sort(c(`Anchor Name`, `Library Name`))[1],
           drug2 = sort(c(`Anchor Name`, `Library Name`))[2])







M1<-matrix(rnorm(1000),nrow=2)
M2<-matrix(rnorm(100),nrow=2)
ee = M1
y = ee[1,] + 10
x = ee[2,] + 10

plot(x, y)
plot(scale(x) - abs(scale(y)), y)
