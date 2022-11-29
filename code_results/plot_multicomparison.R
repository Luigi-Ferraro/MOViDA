########################################################################################################################
############################## Libraries
########################################################################################################################
library(ggplot2)
library(readr)
library(tidyverse)
library(caret)
library(scales)
library(ggridges)
theme_set(theme_minimal())
########################################################################################################################
############################## Inputs
########################################################################################################################
save_bool = F
max_value = 1.2

# lst1 = list(model_directory = "model1_falseweighted_80eps_300e_mmse_newdrfeat",
#             data_directory = "data_newgenes_ge_ont_drfeat",
#             exp_directory = "all_tests_newgenes_newdrfeat",
#             name_displayed = "MOViDA",
#             max_value = 1.2,
#             data_test_file = "drugcell_test")
#
# lst2 = list(model_directory = "model1_falseweighted_80eps_300e_mmse_newdr",
#             data_directory = "data_newgenes_ge_ont_drfeat",
#             exp_directory = "all_tests_newgenes_newdrfeat",
#             name_displayed = "MOViDA_wo_drfeat",
#             max_value = 1.2,
#             data_test_file = "drugcell_test")

lst1 = list(model_directory = "model1_falseweighted_80eps_300e_mmse",
            data_directory = "data_newgenes_ge_ont",
            exp_directory = "all_tests_newgenes",
            name_displayed = "MOViDA",
            max_value = 1.2,
            data_test_file = "drugcell_test")

lst2 = list(model_directory = "pretrained",
            data_directory = "data",
            exp_directory = "all_tests_newgenes",
            name_displayed = "DrugCell",
            max_value = NA,
            data_test_file = "drugcell_test_geont")



#lst_inputs = list(MOViDA = lst1, MOViDA_wo_drfeat = lst2, DrugCell = lst3)
lst_inputs = list(MOViDA = lst1, DrugCell = lst2)








lst_inputs = list(
                  newdr = list(model_directory = "movida_newdr",
                             data_directory = "data_newgenes_ge_ont_drfeat",
                             exp_directory = "all_tests_newgenes_newdrfeat",
                             name_displayed = "newdr",
                             max_value = 1.2,
                             data_test_file = "drugcell_test"),
                  newdr_pcfp = list(model_directory = "movida_newdr_pcfp",
                             data_directory = "data_newgenes_ge_ont_drfeat",
                             exp_directory = "all_tests_newgenes_newdrfeat",
                             name_displayed = "newdr_pcfp",
                             max_value = 1.2,
                             data_test_file = "drugcell_test"),
                  newdrfeat_all = list(model_directory = "movida_newdrfeat_all",
                             data_directory = "data_newgenes_ge_ont_drfeat",
                             exp_directory = "all_tests_newgenes_newdrfeat",
                             name_displayed = "newdrfeat_all",
                             max_value = 1.2,
                             data_test_file = "drugcell_test"),
                  newdrfeat_all_pcfp = list(model_directory = "movida_newdrfeat_all_pcfp",
                             data_directory = "data_newgenes_ge_ont_drfeat",
                             exp_directory = "all_tests_newgenes_newdrfeat",
                             name_displayed = "newdrfeat_all_pcfp",
                             max_value = 1.2,
                             data_test_file = "drugcell_test"),
                  newdrfeat_all_wo_MF = list(model_directory = "movida_newdrfeat_all_wo_MF",
                                              data_directory = "data_newgenes_ge_ont_drfeat",
                                              exp_directory = "all_tests_newgenes_newdrfeat",
                                              name_displayed = "newdrfeat_all_wo_MF",
                                              max_value = 1.2,
                                              data_test_file = "drugcell_test"),
                  DrugCell = list(model_directory = "pretrained",
                             data_directory = "data",
                             exp_directory = "all_tests_newgenes_newdrfeat",
                             name_displayed = "DrugCell",
                             max_value = NA,
                             data_test_file = "drugcell_test")
                  # all_16k = list(model_directory = "model1_falseweighted_80eps_300e_mmse_newdrfeat_all",
                  #                data_directory = "data_newgenes_ge_ont_drfeat",
                  #                exp_directory = "all_tests_newgenes_newdrfeat",
                  #                name_displayed = "all_16k",
                  #                max_value = 1.2,
                  #                data_test_file = "drugcell_test")
                  )




























#
# lst_inputs = list(cv1 = list(model_directory = "model1_falseweighted_80eps_300e_mmse_cv1",
#                              data_directory = "cv1",
#                              exp_directory = "all_tests_newgenes_newdrfeat",
#                              name_displayed = "cv1",
#                              max_value = 1.2,
#                              data_test_file = "drugcell_test"),
#                   cv2 = list(model_directory = "model1_falseweighted_80eps_300e_mmse_cv2",
#                              data_directory = "cv2",
#                              exp_directory = "all_tests_newgenes_newdrfeat",
#                              name_displayed = "cv2",
#                              max_value = 1.2,
#                              data_test_file = "drugcell_test"),
#                   cv3 = list(model_directory = "model1_falseweighted_80eps_300e_mmse_cv3",
#                              data_directory = "cv3",
#                              exp_directory = "all_tests_newgenes_newdrfeat",
#                              name_displayed = "cv3",
#                              max_value = 1.2,
#                              data_test_file = "drugcell_test"),
#                   cv4 = list(model_directory = "model1_falseweighted_80eps_300e_mmse_cv4",
#                              data_directory = "cv4",
#                              exp_directory = "all_tests_newgenes_newdrfeat",
#                              name_displayed = "cv4",
#                              max_value = 1.2,
#                              data_test_file = "drugcell_test"),
#                   cv5 = list(model_directory = "model1_falseweighted_80eps_300e_mmse_cv5",
#                              data_directory = "cv5",
#                              exp_directory = "all_tests_newgenes_newdrfeat",
#                              name_displayed = "cv5",
#                              max_value = 1.2,
#                              data_test_file = "drugcell_test"))
#












########################################################################################################################
############################## Create model function
############################## model_directory, data_directory, exp_directory, name_displayed,
############################## max_value, data_test_file
########################################################################################################################

create_model_df = function(model_directory, data_directory, exp_directory,
                           name_displayed, max_value = NA, data_test_file = "drugcell_test") {
    print(name_displayed)
    model_test_data <- read_delim(paste0(data_directory, "/", data_test_file,".txt"),
                                  "\t", escape_double = FALSE, col_names = FALSE,
                                  trim_ws = TRUE, col_types = cols())
    model_all_data <- read_delim(paste0(data_directory, "/drugcell_all.txt"),
                                 "\t", escape_double = FALSE, col_names = FALSE,
                                 trim_ws = TRUE, col_types = cols())
    model_test_pred <- read_csv(paste0("launch/", exp_directory, "/", model_directory, "/Result/drugcell.predict"),
                                col_names = FALSE, col_types = cols())
    model_all_pred <- read_csv(paste0("launch/", exp_directory, "/", model_directory, "/Result2/drugcell.predict"),
                               col_names = FALSE, col_types = cols())

    if (!is.na(max_value)) {
        model_test_data$X3[model_test_data$X3 >= max_value] <- max_value - 0.000001
        model_test_pred$X1[model_test_pred$X1 < 0] <- 0
        model_test_pred$X1[model_test_pred$X1 >= max_value] <- max_value - 0.000001
        model_all_data$X3[model_all_data$X3 >= max_value] <- max_value - 0.000001
        model_all_pred$X1[model_all_pred$X1 < 0] <- 0
        model_all_pred$X1[model_all_pred$X1 >= max_value] <- max_value - 0.000001
    }

    model_df_test <- data.frame(AUC = model_test_data$X3, AUC_pred = model_test_pred$X1, class = floor(model_test_data$X3 * 10),
                                class_pred = floor(model_test_pred$X1 * 10), ccl = model_test_data$X1, cmp = model_test_data$X2, stringsAsFactors = TRUE)
    model_df_test$type <- name_displayed

    model_df_all <- data.frame(AUC = model_all_data$X3, AUC_pred = model_all_pred$X1, class = floor(model_all_data$X3 * 10),
                               class_pred = floor(model_all_pred$X1 * 10), ccl = model_all_data$X1, cmp = model_all_data$X2, stringsAsFactors = TRUE)
    model_df_all$type <- name_displayed

    model_df_test[model_df_test$class >= 10, "class"] <- Inf
    model_df_test[model_df_test$class_pred >= 10, "class_pred"] <- Inf
    model_df_all[model_df_all$class >= 10, "class"] <- Inf
    model_df_all[model_df_all$class_pred >= 10, "class_pred"] <- Inf

    list(model_df_test = model_df_test, model_df_all = model_df_all)
}


########################################################################################################################
############################## Create model_df for each model
########################################################################################################################
lst_models_df = lapply(lst_inputs, function(x) do.call(create_model_df, x))




########################################################################################################################
############################## Plot distribution samples in classes before and after thresholding
########################################################################################################################

all_data_count <- read_delim(paste0(lst1$data_directory, "/drugcell_all.txt"),
                             "\t", escape_double = FALSE, col_names = FALSE,
                             trim_ws = TRUE, col_types = cols())
all_data_count$class <- floor(all_data_count$X3 * 10)
all_data_count_tmp <- all_data_count %>% group_by(class) %>% count()
p = ggplot(all_data_count_tmp, aes(x=class, y=n)) +
    geom_bar(stat = "identity", fill = 'skyblue', color = 'grey30') +
    theme_bw() +
    labs(title = 'Histogram of scores in ALL', x = 'Scores', y = "counts log10 scale") +
    geom_text(aes(label=n), vjust=-0.3) +
    scale_y_log10()

if (save_bool) {
    ggsave(
        "BarPlot_DataUnbalance_before_tr.pdf",
        plot = last_plot(),
        device = "pdf",
        path = "./code_results/comparisons/",
        scale = 1,
        width = NA,
        height = NA
    )
} else {
    print(p)
}



all_data_count$X3[all_data_count$X3 >= max_value] <- max_value - 0.000001
all_data_count$class <- floor(all_data_count$X3 * 10)
all_data_count_tmp <- all_data_count %>% group_by(class) %>% count()
p = ggplot(all_data_count_tmp, aes(x=class, y=n)) +
    geom_bar(stat = "identity", fill = 'skyblue', color = 'grey30') +
    theme_bw() +
    labs(title = 'Histogram of scores in ALL', x = 'Scores', y = "counts log10 scale") +
    geom_text(aes(label=n), vjust=-0.3) +
    scale_y_log10()

if (save_bool) {
    ggsave(
        "BarPlot_DataUnbalance.pdf.pdf",
        plot = last_plot(),
        device = "pdf",
        path = "./code_results/comparisons/",
        scale = 1,
        width = NA,
        height = NA
    )
} else {
    print(p)
}





########################################################################################################################
############################## Box and ridge plots
########################################################################################################################

all_tests = do.call("rbind", lapply(lst_models_df, function(x) x[[1]]))
rownames(all_tests) = NULL



tmp <- all_tests
tmp$class = as.factor(tmp$class)
tmp$class_pred = as.factor(tmp$class_pred)

p = ggplot(tmp, aes(x= class, y= AUC_pred, color = type)) +
    geom_boxplot(outlier.size = 0.5)+
    #xlim(0, 13) +
    ylim(0, 1.2)  +
    labs(title = 'Score classes vs Prediction TEST', x = 'Score classes', y = 'Prediction') +
    theme_bw()


if (save_bool) {
    ggsave(
        "BoxPlot_score_classes.pdf",
        plot = last_plot(),
        device = "pdf",
        path = "./code_results/comparisons/",
        scale = 1,
        width = NA,
        height = NA
    )
} else {
    print(p)
}

p = ggplot(tmp, aes(x= class_pred, y= AUC, color = type)) +
    geom_boxplot(outlier.size = 0.5)+
    #xlim(0, 1.5) +
    ylim(0, 1.2)  +
    labs(title = 'Prediction classes vs Score TEST', x = 'Prediction classes', y = 'Score') +
    theme_bw()

if (save_bool) {
    ggsave(
        "BoxPlot_prediction_classes.pdf",
        plot = last_plot(),
        device = "pdf",
        path = "./code_results/comparisons/",
        scale = 1,
        width = NA,
        height = NA
    )
} else {
    print(p)
}


p = ggplot(tmp) +
    geom_density_ridges(aes(x = AUC_pred, y = class,
                            group = interaction(type, class),
                            fill = type,
                            vline_color = type),
                        vline_size = 1,
                        quantile_lines = TRUE, quantiles = 2,
                        alpha = 0.3, scale = 2, size = 0.3, rel_min_height = 0.01) +
    theme_bw() +
    scale_x_continuous(breaks = seq(0.0, 1.1, by = 0.1), limits = c(0, 1.11))

if (save_bool) {
    ggsave(
        "Density_ridges_score_classes.pdf",
        plot = last_plot(),
        device = "pdf",
        path = "./code_results/comparisons/",
        scale = 1,
        width = NA,
        height = NA
    )
} else {
    print(p)
}




########################################################################################################################
############################## MSE per class plots
########################################################################################################################

MSE <- function(x, y) {
    mean((x - y)^2)
}

mse_df <- tmp %>% group_by(class, type) %>% summarise(mse=MSE(AUC_pred, AUC))


p = ggplot(mse_df, aes(x = class, y = mse, color = type, group = type)) +
    geom_line() +
    ylim(0, 0.12) +
    labs(title = 'MSE calculated for each class TEST', x = 'Score Class', y = 'MSE') +
    theme_bw() +
    stat_smooth(method="lm", formula=y~1, se=FALSE, size = 0.3, linetype="dashed")


if (save_bool) {
    ggsave(
        "MSE_score_class_comparison.pdf",
        plot = last_plot(),
        device = "pdf",
        path = "./code_results/comparisons/",
        scale = 1,
        width = NA,
        height = NA
    )
} else {
    print(p)
}


mse_df <- tmp %>% group_by(class_pred, type) %>% summarise(mse=MSE(AUC_pred, AUC))

p = ggplot(mse_df, aes(x = class_pred, y = mse, color = type, group = type)) +
    geom_line() +
    ylim(0, 0.15) +
    labs(title = 'MSE calculated for each prediction class TEST', x = 'Prediction Class', y = 'MSE') +
    theme_bw()


if (save_bool) {
    ggsave(
        "MSE_prediction_class_comparison.pdf",
        plot = last_plot(),
        device = "pdf",
        path = "./code_results/comparisons/",
        scale = 1,
        width = NA,
        height = NA
    )
} else {
    print(p)
}



########################################################################################################################
############################## Confusion Matrices
########################################################################################################################


plot_confusion_matrix <- function(cm_df, save_bool = F, color_high = "red", color_low = "lightblue1") {
    name_displ = unique(cm_df$type)
    cm_df$class = as.factor(cm_df$class)
    cm_df$class_pred = as.factor(cm_df$class_pred)
    cm <- confusionMatrix(data = cm_df$class_pred, reference = cm_df$class)
    tab <- as.data.frame(cm$table)
    freq_tab <- tab %>% group_by(Reference) %>% mutate(percent = Freq/sum(Freq))
    p <- ggplot(freq_tab, aes(Reference, Prediction)) +
        geom_tile(aes(fill = percent)) +
        geom_text(aes(label = Freq)) + #paste(round(percent * 100, 2), " %")))  +
        labs(title = paste('Confusion Matrix ALL', name_displ), x = 'Score Class', y = 'Prediction Class') +
        theme_bw() +
        scale_fill_gradientn(
            colors=c(color_low, color_high),
            values=rescale(c(0,0.7)),
            limits=c(0,0.7)
        )


    if (save_bool) {
        ggsave(
            paste0("CM_", name_displ, ".pdf"),
            plot = last_plot(),
            device = "pdf",
            path = "./code_results/comparisons/",
            scale = 1,
            width = NA,
            height = NA
        )
    } else {
        print(p)
    }
}



lapply(lst_models_df, function(x) plot_confusion_matrix(x[[1]], save_bool))




########################################################################################################################
############################## Tables
########################################################################################################################

lst_names = as.vector(sapply(lst_inputs, function(x) x$name_displayed))
res_tab = lapply(lst_models_df, function(x)
                            x[[1]] %>% group_by(class) %>%
                            summarise(mse = mean( (AUC - AUC_pred)**2 )))
for (i in 1:length(lst_names)) {
    colnames(res_tab[[i]])[2] = paste0("mse_", lst_names[i])
}

mse_byclass_all = Reduce(function(x, y) merge(x, y, by = "class"), res_tab)
mse_byclass_all$mse_mean = apply(mse_byclass_all[, 2:(length(lst_names)+1)], 1, mean)
mse_byclass_all



res_tab = lapply(lst_models_df, function(x) x[[1]] %>% group_by(class) %>%
                     summarise(mse = mean( (AUC - AUC_pred)**2 )) %>%
                     summarise(mmse = mean(mse)))

res_bind1 = as.data.frame(do.call(rbind.data.frame, res_tab))
rownames(res_bind1) = names(res_tab)


res_tab = lapply(lst_models_df, function(x) x[[1]] %>%
                     summarise(cor = cor(AUC, AUC_pred),
                               cor_t = cor.test(AUC, AUC_pred)$statistic,
                               cor_pval = cor.test(AUC, AUC_pred)$p.value,
                               scor = cor(AUC, AUC_pred, method = "spearman"),
                               scor_t = cor.test(AUC, AUC_pred, method = "spearman")$statistic,
                               scor_pval = cor.test(AUC, AUC_pred, method = "spearman")$p.value,
                               cor_class = cor(class, class_pred),
                               scor_class = cor(class, class_pred, method = "spearman")))

res_bind2 = as.data.frame(do.call(rbind.data.frame, res_tab))
rownames(res_bind2) = names(res_tab)

round(cbind(res_bind1, res_bind2), digits = 2)

noninfite_op = function(lst_models_df) {
    for ( i in 1:length(lst_models_df) ) {
        tmpc = lst_models_df[[i]][[1]]$class
        tmpc[is.infinite(tmpc)] = 10
        tmpcp = lst_models_df[[i]][[1]]$class_pred
        tmpcp[is.infinite(tmpcp)] = 10
        message(names(lst_models_df)[i])
        print(paste("Pearson: ", cor(tmpc, tmpcp)))
        print(paste("Spearman: ", cor(tmpc, tmpcp, method = "spearman")))

    }
}
noninfite_op(lst_models_df)



n_class_0 = sum(lst_models_df[[1]][[1]]$class == 0)


n_class_0 = 100
res_sampl = c()
for (x in 1:1000) {
    tmp = lapply(lst_models_df, function(x) {tmp = x[[1]] %>%
        group_by(class) %>%
        sample_n(n_class_0) %>%
        ungroup() %>%
        summarise(cor = cor(AUC, AUC_pred),
                  scor = cor(AUC, AUC_pred, method = "spearman"))})
    res_sampl[[x]] = tmp
}

lst_res = list()
for (i in names(lst_models_df)) {
    lst_res[[i]] = data.frame(matrix(0, 0, 2))
    colnames(lst_res[[i]]) = c("cor", "scor")
    for (x in 1:length(res_sampl)) {
        lst_res[[i]] = rbind(lst_res[[i]], res_sampl[[x]][[i]])
    }
}


lapply(lst_res, function(x) round(apply(x, 2, mean), digits = 2))

