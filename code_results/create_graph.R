library(igraph)
library(networkD3)
library(visNetwork)

#HAEMATOPOIETIC_AND_LYMPHOID_TISSUE
drug_sel = "Bleomycin sulfate"
family_sel = "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"
head_x = 5

ont_tree <- read_delim("data_newgenes_ge_ont/ont_tree.txt",
                       "\t", escape_double = FALSE, col_names = FALSE,
                       trim_ws = TRUE)
colnames(ont_tree) = c("father", "child")

if (family_sel == "") {
    sel_dr_fam = important_go_gene_in_drugfamily_all %>% subset(drug == drug_sel)
} else {
    sel_dr_fam = important_go_gene_in_drugfamily %>% subset(drug == drug_sel & family == family_sel)
}
sel_dr_fam = sel_dr_fam %>% head(head_x)


tmp = important_go_gene_in_drugfamily %>% group_by(drug, family) %>%
    filter(row_number() %in% 1:head_x)


gos_check = sel_dr_fam$go#union(ont_tree$father, ont_tree$child)

visited = c()
df_edges_go = data.frame()
nodes_layers = data.frame(label = gos_check, level = 5)
while (length(gos_check) > 0) {
    go = gos_check[1]
    gos_check = gos_check[-1]
    if (!go %in% visited) {
        visited = c(visited, go)
        go_fathers = ont_tree %>% subset(child == go)
        gos_check = c(go_fathers$father, gos_check)
        df_edges_go = rbind(df_edges_go, go_fathers)
        go_fathers$level = nodes_layers[nodes_layers$label == go, "level"] - 1
        go_fathers = go_fathers[, c(1,3)]
        colnames(go_fathers) = c("label", "level")
        go_fathers = go_fathers %>% subset(!label %in% nodes_layers$label)
        nodes_layers = rbind(nodes_layers, go_fathers)
    }
}
nodes_layers[nodes_layers$label == "GO:0008150", "level"] = 1#max(nodes_layers$level) + 1
nodes_layers$color = "lightblue"

sel_dr_fam_genes = sel_dr_fam %>% subset(   genes_mut != "" |
                                            genes_amp != "" |
                                            genes_del != "")
for (i in 1:nrow(sel_dr_fam_genes)) {
    row = sel_dr_fam_genes[i,]
    if (row$genes_mut != "") {
        genes = as.vector(str_split(row$genes_mut, ",", simplify = T))
        df_edges_go = rbind(df_edges_go, data.frame(father = row$go, child = genes))
        nodes_layers = rbind(nodes_layers,
                             data.frame(label = genes, level = 100, color = "red"))
    }
    if (row$genes_amp != "") {
        genes = as.vector(str_split(row$genes_amp, ",", simplify = T))
        df_edges_go = rbind(df_edges_go, data.frame(father = row$go, child = genes))
        nodes_layers = rbind(nodes_layers,
                             data.frame(label = genes, level = 100, color = "green"))
    }
    if (row$genes_del != "") {
        genes = as.vector(str_split(row$genes_del, ",", simplify = T))
        df_edges_go = rbind(df_edges_go, data.frame(father = row$go, child = genes))
        nodes_layers = rbind(nodes_layers,
                             data.frame(label = genes, level = 100, color = "purple"))
    }
}

nodes_layers = unique(nodes_layers)

nodes_layers$id = 1:nrow(nodes_layers)
backward = c("GO:0008150")
i = 1
while (length(backward) > 0) {
    backward = df_edges_go %>% subset(father %in% backward) %>% arrange(child)
    backward = backward$child
    i = i + 1
    nodes_layers[nodes_layers$label %in% backward, "level"] = i
}

colnames(nodes_layers)[1] = "father"
df_edges = df_edges_go %>% left_join(nodes_layers, by = "father")
df_edges = df_edges[, c(5,2)]
colnames(df_edges)[1] = "father"

colnames(nodes_layers)[1] = "child"
df_edges = df_edges %>% left_join(nodes_layers, by = "child")
df_edges = df_edges[, c(1,5)]
colnames(df_edges)[2] = "child"

colnames(nodes_layers)[1] = "label"
colnames(df_edges) = c("from", "to")

df_edges = df_edges %>% arrange(to)
nodes_layers = nodes_layers %>% arrange(label)

gos_descr_space = within(gos_descr, descr_space <- gsub("(([-/A-Za-z1-9.,']+\\s){3})","\\1\n ", descr))

nodes_layers_descr = nodes_layers
colnames(nodes_layers_descr)[1] = "go"
nodes_layers_descr = nodes_layers_descr %>% left_join(gos_descr_space, by = "go")
nodes_layers_descr$descr_space <- ifelse(is.na(nodes_layers_descr$descr_space),
                               nodes_layers_descr$go, nodes_layers_descr$descr_space)
nodes_layers_descr = nodes_layers_descr[, c(6,2,4,3)]
colnames(nodes_layers_descr)[1] = "label"


visNetwork(nodes_layers_descr, df_edges) %>%
    visEdges(smooth = T) %>%
    visNodes() %>%
    visHierarchicalLayout(direction = "LR", levelSeparation = 1000)





