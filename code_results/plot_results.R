
library(ggplot2)
library(rlang)
library(dplyr)
library(magrittr)
library(MLmetrics)
library(readr)
library(caret)
library(pROC)

dirname <- "model1_20eps_100e"
dirdataname <- "data_newgenes_ge_ont"
bool_max_value = TRUE
bool_train_info = TRUE
max_value <- 1.2
direxp <- "all_tests_newgenes"

save_bool <- FALSE


calc_f1_score <- function(df) {
    cm <- confusionMatrix(data = df$class_pred, reference = df$class)
    cm_tab <- cm$table
    precision <- diag(cm_tab) / rowSums(cm_tab)
    recall <- diag(cm_tab) / colSums(cm_tab)
    f1_pclass <- 2 * (precision * recall)/(precision + recall)
    f1_pclass
}

plot_res <- function(dirname, dirdataname, direxp, save_bool, max_value = 1.2, bool_max_value = TRUE, bool_train_info = TRUE) {
    print(dirname)
    print(max_value)
    if (bool_train_info) {
        train_info <- read_csv(paste0("./launch/", direxp, "/", dirname, "/train_info.csv"), col_types = cols())
    }
    drugcell2 <- read_csv(paste0("./launch/", direxp, "/", dirname, "/Result2/drugcell.predict"),
                          col_names = FALSE, col_types = cols())
    drugcell <- read_csv(paste0("./launch/", direxp, "/", dirname, "/Result/drugcell.predict"),
                         col_names = FALSE, col_types = cols())
    drugcell_all <- read_delim(paste0(dirdataname, "/drugcell_all.txt"),
                               "\t", escape_double = FALSE, col_names = FALSE,
                               trim_ws = TRUE, col_types = cols())
    drugcell_test <- read_delim(paste0(dirdataname, "/drugcell_test.txt"),
                                "\t", escape_double = FALSE, col_names = FALSE,
                                trim_ws = TRUE, col_types = cols())


    #h <- hist(drugcell_all$X3, breaks = 25)
    #print(h)
    if (bool_max_value) {
        drugcell_all$X3[drugcell_all$X3 >= max_value] <- max_value - 0.000001
        drugcell_test$X3[drugcell_test$X3 >= max_value] <- max_value - 0.000001
        drugcell$X1[drugcell$X1 < 0] <- 0
        drugcell2$X1[drugcell2$X1 < 0] <- 0
        drugcell$X1[drugcell$X1 >= max_value] <- max_value - 0.000001
        drugcell2$X1[drugcell2$X1 >= max_value] <- max_value - 0.000001
    }
    #h <- hist(drugcell_all$X3, breaks = max_value)
    #print(h)
    #h <- hist(drugcell_all$X3, breaks = 10)


    test <- drugcell_test
    df = data.frame(AUC = test$X3, AUC_pred = drugcell$X1, class = floor(test$X3 * 10), class_pred = floor(drugcell$X1 * 10), cmp = test$X2, stringsAsFactors = TRUE)
    df$class = as.factor(df$class)
    df$class_pred = as.factor(df$class_pred)

    all <- drugcell_all
    df_all = data.frame(AUC = all$X3, AUC_pred = drugcell2$X1, class = floor(all$X3 * 10), class_pred = floor(drugcell2$X1 * 10), cmp = all$X2, stringsAsFactors = TRUE)
    df_all$class = as.factor(df_all$class)
    df_all$class_pred = as.factor(df_all$class_pred)


    if (save_bool) pdf(paste0("./launch/", direxp, "/", dirname, "/results_plots.pdf"))
    #plot(df$AUC, df$AUC_pred, cex = 1.1)

    #cor(df$AUC, df$AUC_pred)
    #abline(0,1, col = "red")
    df_cap <- df
    df_cap$class <- as.integer(df_cap$class) - 1
    df_cap$class_pred <- as.integer(df_cap$class_pred) - 1
    df_cap[df_cap$class >= 10, "class"] <- Inf
    df_cap[df_cap$class_pred >= 10, "class_pred"] <- Inf
    df_cap$class <- as.factor(df_cap$class)
    df_cap$class_pred <- as.factor(df_cap$class_pred)
    df <- df_cap

    df_cap <- df_all
    df_cap$class <- as.integer(df_cap$class) - 1
    df_cap$class_pred <- as.integer(df_cap$class_pred) - 1
    df_cap[df_cap$class >= 10, "class"] <- Inf
    df_cap[df_cap$class_pred >= 10, "class_pred"] <- Inf
    df_cap$class <- as.factor(df_cap$class)
    df_cap$class_pred <- as.factor(df_cap$class_pred)
    df_all <- df_cap

    p = ggplot(df, aes(x = AUC, y = AUC_pred)) +
        geom_point(size = 0.5, alpha = 1/10, color = "blue") +
        xlim(0, 1.25) +
        ylim(0, 1.25) +
        geom_abline() +
        labs(title = 'Score vs Prediction TEST', x = 'Score', y = 'Prediction')


    print(p)

    p = ggplot(df, aes(x = AUC, y = AUC_pred)) +
        geom_bin2d(bins = 500) +
        xlim(0, 1.25) +
        ylim(0, 1.25) +
        scale_fill_continuous(type = "viridis") +
        geom_abline()+
        labs(title = 'Score vs Prediction TEST', x = 'Score', y = 'Prediction')


    print(p)

    p = ggplot(df_all, aes(x = AUC, y = AUC_pred)) +
        geom_bin2d(bins = 500) +
        xlim(0, 1.25) +
        ylim(0, 1.25) +
        scale_fill_continuous(type = "viridis") +
        geom_abline()+
        labs(title = 'Score vs Prediction ALL dataset', x = 'Score', y = 'Prediction')


    print(p)

    df$type <- "test"
    df_all$type <- "all"
    tmp1 <- subset(df, select = c("class", "AUC_pred", "type"))
    tmp2 <- subset(df_all, select = c("class", "AUC_pred", "type"))
    tmp <- rbind(tmp1, tmp2)

    p = ggplot(tmp, aes(x= class, y= AUC_pred, color = type)) +
        geom_boxplot()+
        #xlim(0, 1.5) +
        ylim(0, 1.25) +
        labs(title = 'Score classes vs Prediction', x = 'Score classes', y = 'Prediction')

    print(p)

    tmp2 <- tmp
    tmp2$class <- as.integer(tmp2$class)
    tmp2 <- tmp2[tmp2$class < 10,]
    tmp2$class <- as.factor(tmp2$class)

    p = ggplot(tmp2, aes(x= class, y= AUC_pred, color = type)) +
        geom_violin()+
        #xlim(0, 1.5) +
        ylim(0, 1.5)  +
        labs(title = 'Score classes vs Prediction', x = 'Score classes', y = 'Prediction')

    print(p)

    tmp1 <- subset(df, select = c("class_pred", "AUC", "type"))
    tmp2 <- subset(df_all, select = c("class_pred", "AUC", "type"))
    tmp <- rbind(tmp1, tmp2)

    p = ggplot(tmp, aes(x= class_pred, y= AUC, color = type)) +
        geom_boxplot()+
        #xlim(0, 1.5) +
        ylim(0, 1.25) +
        labs(title = 'Score vs Prediction classes', x = 'Prediction classes', y = 'Score')

    print(p)



    tmp2 <- tmp
    tmp2$class <- as.integer(tmp2$class)
    tmp2 <- tmp2[tmp2$class < 10,]
    tmp2$class <- as.factor(tmp2$class)

    p = ggplot(tmp2, aes(x= class_pred, y= AUC, color = type)) +
        geom_violin()+
        #xlim(0, 1.5) +
        ylim(0, 1.5)  +
        labs(title = 'Score vs Prediction classes', x = 'Prediction classes', y = 'Score')

    print(p)


    p <- ggplot() +
        geom_density(aes(x = df$AUC, colour="score", fill="score"), alpha=.2) +
        geom_vline(data=df, aes(xintercept=mean(AUC), colour="score",linetype="mean"), size=0.5) +
        geom_vline(data=df, aes(xintercept=median(AUC), colour="score",linetype="median"), size=0.5) +
        geom_density(aes(x = df$AUC_pred, colour="pred", fill="pred"), alpha=.2) +
        geom_vline(data=df, aes(xintercept=mean(AUC_pred), colour="pred",linetype="mean"), size=0.5) +
        geom_vline(data=df, aes(xintercept=median(AUC_pred), colour="pred",linetype="median"), size=0.5) +
        labs(title = 'Density of Scores and Predicion on TEST', x = 'Scores', y = 'Density') +
        ylim(0, 10)

    print(p)



    tmp <- df  %>%
        group_by(class) %>%
        summarize(mean = mean(AUC_pred),
                  median = median(AUC_pred))

    p <- ggplot() +
        geom_density(aes(x = df$AUC_pred, colour=df$class, fill=df$class), alpha=.2) +
        geom_vline(data=tmp, aes(xintercept=mean, colour=class, linetype="mean"), size=0.5) +
        geom_vline(data=tmp, aes(xintercept=median, colour=class, linetype="median"), size=0.5) +
        labs(title = 'Density of Scores and Predicion on TEST by classes', x = 'Scores', y = 'Density') +
        ylim(0, 10)

    print(p)
    #ggplot(df, aes(x= cmp, y= AUC_pred)) +
    #    geom_boxplot()


    #grouped <- df %>% group_by(cmp)


    #df_corr = by(df, df$cmp, function(x) cor(x$AUC,x$AUC_pred))
    #hist(df_corr)


    if (bool_train_info) {
        p = ggplot(train_info) +
            geom_line(aes(x= epoch, y= train_corr, color = "train")) +
            geom_line(aes(x= epoch, y= val_corr, color = "val")) +
            ylim(0, 1) +
            labs(title = 'Correlation per epoch on Training and Validation sets', x = 'Epoch', y = 'Pearson correlation')

        print(p)
    }




    if (bool_train_info) {
        p = ggplot(train_info) +
            geom_line(aes(x= epoch, y= train_mse, color = "train")) +
            geom_line(aes(x= epoch, y= val_mse, color = "val")) +
            ylim(0, 1) +
            labs(title = 'MSE per epoch on Training and Validation sets', x = 'Epoch', y = 'MSE')

        print(p)
    }




    mse_groups <- df %>% group_by(class) %>% summarise(mse=MSE(AUC_pred, AUC),
                                                       pearson = cor(AUC, AUC_pred),
                                                       spearman = cor(AUC, AUC_pred, method = "spearman"))
    #print(unique(mse_groups$class))
    mse_groups$class <- as.integer(mse_groups$class) - 1
    mse_groups$class <- mse_groups$class
    if (bool_max_value) {
        mse_groups_lim <- mse_groups[mse_groups$class < max_value * 10, ]
    } else {
        mse_groups_lim <- mse_groups
    }
    mse_groups_lim$class = as.factor(mse_groups_lim$class)
    mse_groups_lim$f1 <- calc_f1_score(df)
    #print(unique(mse_groups_lim$class))

    mse_groups2 <- df_all %>% group_by(class) %>% summarise(mse=MSE(AUC_pred, AUC),
                                                            pearson = cor(AUC, AUC_pred),
                                                            spearman = cor(AUC, AUC_pred, method = "spearman"))
    #print(unique(mse_groups$class))
    mse_groups2$class <- as.integer(mse_groups2$class) - 1
    mse_groups2$class <- mse_groups2$class
    if (bool_max_value) {
        mse_groups_lim2 <- mse_groups2[mse_groups2$class < max_value * 10, ]
    } else {
        mse_groups_lim2 <- mse_groups2
    }
    mse_groups_lim2$class = as.factor(mse_groups_lim2$class)
    mse_groups_lim2$f1 <- calc_f1_score(df_all)
    #print(unique(mse_groups_lim$class))

    mse_groups_lim$type <- "test"
    mse_groups_lim2$type <- "all"
    tmp <- rbind(mse_groups_lim, mse_groups_lim2)

    p = ggplot(tmp, aes(x = class, y = mse, color = type, group = type)) +
        geom_line() +
        ylim(0, 0.15) +
        labs(title = 'MSE calculated for each class', x = 'Score Class', y = 'MSE')

    print(p)

    p = ggplot(tmp, aes(x = class, y = pearson, color = type, group = type)) +
        geom_line() +
        ylim(0, 0.5) +
        labs(title = 'Pearson correlation calculated for each class', x = 'Score Class', y = 'Pearson')

    print(p)

    p = ggplot(tmp, aes(x = class, y = spearman, color = type, group = type)) +
        geom_line() +
        ylim(0, 0.5) +
        labs(title = 'Spearman correlation calculated for each class', x = 'Score Class', y = 'Spearman')

    print(p)



    p = ggplot(tmp, aes(x = class, y = f1, color = type, group = type)) +
        geom_line() +
        ylim(0, 1) +
        labs(title = 'F1 score calculated for each class', x = 'Score Class', y = 'F1 score')

    print(p)





    mse_groups <- df %>% group_by(class_pred) %>% summarise(mse=MSE(AUC_pred, AUC))
    #print(unique(mse_groups$class))
    mse_groups$class_pred <- as.integer(mse_groups$class_pred)
    mse_groups$class_pred <- mse_groups$class_pred - 1
    if (bool_max_value) {
        mse_groups_lim <- mse_groups[mse_groups$class_pred < max_value * 10, ]
    } else {
        mse_groups_lim <- mse_groups
    }
    mse_groups_lim$class_pred = as.factor(mse_groups_lim$class_pred)
    #print(unique(mse_groups_lim$class))

    mse_groups2 <- df_all %>% group_by(class_pred) %>% summarise(mse=MSE(AUC_pred, AUC))
    #print(unique(mse_groups$class))
    mse_groups2$class_pred <- as.integer(mse_groups2$class_pred)
    mse_groups2$class_pred <- mse_groups2$class_pred - 1
    if (bool_max_value) {
        mse_groups_lim2 <- mse_groups2[mse_groups2$class_pred < max_value * 10, ]
    } else {
        mse_groups_lim2 <- mse_groups2
    }
    mse_groups_lim2$class_pred = as.factor(mse_groups_lim2$class_pred)
    #print(unique(mse_groups_lim$class))

    mse_groups_lim$type <- "test"
    mse_groups_lim2$type <- "all"
    tmp <- rbind(mse_groups_lim, mse_groups_lim2)

    p = ggplot(tmp, aes(x = class_pred, y = mse, color = type, group = type)) +
        geom_line() +
        ylim(0, 0.15) +
        labs(title = 'MSE calculated for each Prediction class', x = 'Prediction Class', y = 'MSE')

    print(p)



    df_cap <- df
    df_cap$class <- as.integer(df_cap$class) - 1
    df_cap$class_pred <- as.integer(df_cap$class_pred) - 1
    df_cap[df_cap$class >= 10, "class"] <- Inf
    df_cap[df_cap$class_pred >= 10, "class_pred"] <- Inf
    df_cap$class <- as.factor(df_cap$class)
    df_cap$class_pred <- as.factor(df_cap$class_pred)
    cm <- confusionMatrix(data = df_cap$class_pred, reference = df_cap$class)
    tab <- as.data.frame(cm$table)
    freq_tab <- tab %>% group_by(Reference) %>% mutate(percent = Freq/sum(Freq))
    p <- ggplot(freq_tab, aes(Reference, Prediction)) +
        geom_tile(aes(fill = percent)) +
        geom_text(aes(label = Freq)) + #paste(round(percent * 100, 2), " %")))  +
        labs(title = 'Confusion Matrix', x = 'Score Class', y = 'Prediction Class')

    print(p)

    p <- ggplot(df, aes(sample=AUC_pred - AUC, color = class), alpha = 0.2) +
        stat_qq() +
        stat_qq_line() +
        #xlim(-1, 1) +
        #ylim(-0.25, 0.5) +
        labs(title = 'QQ-plot')

    print(p)

    p <- ggplot() +
        geom_density(aes(x = df$AUC_pred - df$AUC, colour = df$class, fill = df$class), alpha=.1) +
        xlim(-1, 1) +
        ylim(0, 15) +
        labs(title = 'Density of Residuals', x = 'Residuals')

    print(p)

    if (save_bool) dev.off()

    print("OK")


    #############################################################################################################





    df_app <- df
    df_app$class <- as.numeric(df_app$class)
    df_app$class_pred <- as.numeric(df_app$class_pred)
    df_app$class[df_app$class <= 2] <- 0
    df_app$class[df_app$class > 2] <- 1
    df_app$class_pred[df_app$class_pred <= 2] <- 0
    df_app$class_pred[df_app$class_pred > 2] <- 1



    print(paste("F1 score bclass", round(F1_Score(y_pred = df_app$class_pred, y_true = df_app$class), 4)))
    #print(paste("F1 score", round(F1_Score(y_pred = df$class_pred, y_true = df$class),4)))
    rss <- sum((df$AUC_pred - df$AUC) ^ 2)
    tss <- sum((df$AUC - mean(df$AUC)) ^ 2)
    rsq <- 1 - rss/tss

    print(paste("R-squared", round(rsq,4)))


    print(paste("MSE", round(MSE(df$AUC_pred, df$AUC), 4)))

    print(df %>%
        group_by(class) %>%
        summarize(
                  mse = MSE(AUC, AUC_pred),
                  pearson_cor = cor(AUC, AUC_pred),
                  pearson_cor_pvalue = cor.test(AUC, AUC_pred)$p.value,
                  spearman_cor = cor(AUC, AUC_pred, method = "spearman"),
                  spearman_cor_pvalue = cor.test(AUC, AUC_pred, method = "spearman")$p.value,
                  wilcoxon_less = wilcox.test(AUC, AUC_pred, paired = TRUE, alternative = "less")$p.value,
                  wilcoxon_greater = wilcox.test(AUC, AUC_pred, paired = TRUE, alternative = "greater")$p.value,
                  wilcoxon_twosided = wilcox.test(AUC, AUC_pred, paired = TRUE)$p.value
                  )
        )





    print(paste("Pearson", round(cor(df$AUC, df$AUC_pred), 4)))
    print(paste("Spearman", round(cor(df$AUC, df$AUC_pred, method = "spearman"), 4)))
    print(paste("Wilcoxon less", wilcox.test(df$AUC, df$AUC_pred, paired = TRUE, alternative = "less")$p.value))
    print(paste("Wilcoxon greater", wilcox.test(df$AUC, df$AUC_pred, paired = TRUE, alternative = "greater")$p.value))
    print(paste("Wilcoxon two sided", wilcox.test(df$AUC, df$AUC_pred, paired = TRUE)$p.value))
    #print(ConfusionDF(y_pred = df_app$class_pred, y_true = df_app$class))
    #print(confusionMatrix(data = df$class_pred, reference = df$class))

    #print(paste("F1 score", round(F1_Score(y_pred = df$class_pred, y_true = df$class),4)))
    cm <- confusionMatrix(data = df$class_pred, reference = df$class)
    cm_tab <- cm$table
    precision <- diag(cm_tab) / rowSums(cm_tab)
    recall <- diag(cm_tab) / colSums(cm_tab)
    f1_pclass <- 2 * (precision * recall)/(precision + recall)
    print(round(f1_pclass, 4))
    print(paste("F1-score:", round(mean(f1_pclass), 4)))
    #print(mse_groups_lim)
    cat("\n")

}



couples <- rbind(

    c("model1_falseweighted_80eps_300e_mmse_newdrfeat", "data_newgenes_ge_ont_drfeat", "all_tests_newgenes_newdrfeat", 1.2, TRUE, TRUE)

)

save_bool <- TRUE

for (i in 1:nrow(couples)) {
    try(plot_res(couples[i, 1], couples[i, 2], couples[i, 3], save_bool, as.double(couples[i, 4]), couples[i, 5], couples[i, 6]))
}



















