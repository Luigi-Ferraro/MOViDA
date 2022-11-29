library(readr)
drugcell_all <- read_delim("data_newgenes_ge_ont/drugcell_all.txt",
                           "\t", escape_double = FALSE, col_names = FALSE,
                           trim_ws = TRUE)

n_groups = 5
d = sample(1:nrow(drugcell_all))
m = length(d)/n_groups
s = split(d, ceiling(seq_along(d)/m))

for (i in 1:n_groups) {
    c_point = floor(length(s[[i]])/2)
    test_idx = s[[i]][1:c_point]
    val_idx = s[[i]][(c_point+1):length(s[[i]])]
    train_idx = Reduce(c, s[-i])

    test = drugcell_all[test_idx,]
    val = drugcell_all[val_idx,]
    train = drugcell_all[train_idx,]

    write_delim(test, paste0("./cv", i, "/drugcell_test.txt"), delim = "\t", col_names = F)
    write_delim(train, paste0("./cv", i, "/drugcell_train.txt"), delim = "\t", col_names = F)
    write_delim(val, paste0("./cv", i, "/drugcell_val.txt"), delim = "\t", col_names = F)

}
