
library(rcdk)

# Output
movida_fullpath_out <- "C:\\Users\\Emanuele\\Documents\\Trieste\\ProgettiRicerca\\2022_Napoli\\drug_fingerprint.csv"


# Lettura del file con le smiles (due colonne, separatore spazio)
movida_fullpath_smi <- "C:\\Users\\Emanuele\\Documents\\Trieste\\ProgettiRicerca\\2022_Napoli\\drug2ind.smi"

movida_df <- read.delim2(file=movida_fullpath_smi, sep=" ")
nr_drugs <-  nrow(movida_df)


# I codici del file pdf vanno da 0 a 880
# I codici in fps@bits vanno da 1 a 881
# Quindi per risalire a che cosa corrisponde un bit, bisogna togliere uno
# rispetto al codice che si ha in R
# nr_pubchem_des <- 881
pubchem_header <- paste0("PubchemBit", c(0:880))
pubchem_values <- rep(0, 881)


# Definizione matrice con tutti 0
movida_matrix <- matrix(nrow=nr_drugs, ncol=881, 0)
colnames(movida_matrix) <- pubchem_header
rownames(movida_matrix) <- movida_df$drugid

for (d in 1:nrow(movida_df)) {
  drugid <- movida_df$drugid[d]
  smiles <- movida_df$smiles[d]
  mol <- rcdk::parse.smiles(smiles)[[1]]
  fps <- rcdk::get.fingerprint(mol, type='pubchem')
  
  # Assign (1) values corresponding to pubchem bits
  movida_matrix[d, fps@bits] <- 1
}

write.csv2(
  x=movida_matrix,
  file=movida_fullpath_out
)
