import sys
import torch
import numpy    as np
import pandas   as pd
import networkx as nx
import networkx.algorithms.components.connected as nxacc

import warnings
warnings.filterwarnings("ignore")

class DataMOViDA:

    # TRAIN NEEDED INPUTS
    input_dir           = None
    names_mo            = None
    cell2id_mapping     = None
    drug2id_mapping     = None
    mo_gene2idx         = None
    genes_dim           = None
    drugs_dim           = None
    train_data          = None
    test_data           = None
    val_data            = None
    all_data            = None
    ccl_features        = None
    drug_features       = None
    path_act            = None
    term2id_mapping     = None
    hierarchy_graph     = None
    root                = None

    def __init__(self, input_dir, names_mo,
                drug2id, cell2id, mo_gene2id,
                train_file, val_file, test_file, all_file,
                cell2features, drug2features, pathway_act_file,
                ontology, mo_gene2ontology):
        
        self.input_dir          = input_dir + "/"

        self.names_mo           = names_mo

        # Dictionaries with position mapping
        self.cell2id_mapping    = self.load_mapping(self.input_dir + cell2id)
        self.drug2id_mapping    = self.load_mapping(self.input_dir + drug2id)
        self.mo_gene2idx        = [self.load_mapping(self.input_dir + gene2id) for gene2id in mo_gene2id]

        # Genes dims
        self.genes_dim = [len(x) for x in self.mo_gene2idx]

        # Read train val test files
        ### IT ASSUMES THAT YOU HAVE ALREADY SPLITTED IN TRAIN/TEST/VAL
        self.train_data, self.test_data, self.val_data, self.all_data  = \
                [self.load_data(self.input_dir + tfile, self.cell2id_mapping, self.drug2id_mapping) \
                    for tfile in [train_file, test_file, val_file, all_file]]

        # Load cell and drug features
        self.ccl_features       = [np.genfromtxt(self.input_dir + feat_file, delimiter=',') for feat_file in cell2features]
        self.drug_features      = [np.genfromtxt(self.input_dir + feat_file, delimiter=',') for feat_file in drug2features]


        # Drug dims
        self.drugs_dim          = [x.shape[1] for x in self.drug_features]


        # Load pathway activities
        self.path_act, self.term2id_mapping = self.prepare_path_act(self.input_dir + pathway_act_file, self.cell2id_mapping)
        
        # Load ontology 
        self.hierarchy_graph, self.root = self.load_ontology(self.input_dir + ontology)

        # Load association between terms and genes
        self.mo_gene2term = [self.get_go_gene_map(self.input_dir + ont_gene, self.hierarchy_graph, gene2idx)  \
                                for ont_gene,gene2idx  in zip(mo_gene2ontology, self.mo_gene2idx)]


        
    def load_mapping(self, mapping_file):
        '''
        Read a tab file: val \t key
        Return a dictionary
        '''
        with open(mapping_file, "r") as f:
            mapping = {key: int(val) for val,key in (line.rstrip().split() for line in f)}
        return mapping
    
    def load_data(self, file_name, cell2id, drug2id):
        '''
        Read train/val/test file, a tab file (ccl, cmp, auc)
        Return a tuple with a feature tensor and label tensor
        '''
        feature = []
        label = []

        with open(file_name, 'r') as f:
            lines = [line.strip().split('\t') for line in f]
        feature = [(cell2id[tokens[0]], drug2id[tokens[1]]) for tokens in lines]
        label = [float(tokens[2]) for tokens in lines]
        
        return (torch.Tensor(feature), torch.Tensor(label))
    
    def prepare_path_act(self, pathway_act_file, cell2id_mapping):
        '''
        Read a tab file with Pathway Activities
        '''
        df_pact = pd.read_csv(pathway_act_file, delimiter='\t')
        df_pact["cell_id"] = [cell2id_mapping[x] for x in df_pact["CCLE_Name"]]

        # the index of this dataframe will be the mapping of cell lines
        df_pact = df_pact.set_index("cell_id").drop(columns=["CCLE_Name", "COSMICID"], errors='ignore')

        term2id_mapping = {term : idx for idx,term in enumerate(df_pact.columns)}
        return df_pact.to_numpy(), term2id_mapping

    def load_ontology(self, file_name):
        '''
        Create a graph based on ontology
        Return the graph and the root
        '''
        dG = nx.DiGraph()
        with open(file_name, "r") as f:
            for line in f: 
                line = line.rstrip().split()
                dG.add_edge(line[0], line[1])

        leaves = [n for n in dG.nodes if dG.in_degree(n) == 0]

        uG = dG.to_undirected()
        connected_subG_list = list(nxacc.connected_components(uG))

        if len(leaves) > 1:
            print('There are more than 1 root of ontology. Please use only one root.')
            sys.exit(1)
        if len(connected_subG_list) > 1:
            print('There are more than connected components. Please connect them.')
            sys.exit(1)

        return dG, leaves[0]
    
    def get_go_gene_map(self, ont_file, dG, mapping):
        '''
        Read associations between genes and terms of hierarchy
        Return a dictionary {term : list of associated genes}
        '''
        df = pd.read_csv(ont_file, delimiter="\t",
                        index_col=False, names=["go", "genes"])
        res = {}

        for term in dG.nodes():
            tmp = set(df[df["go"] == term]["genes"])
            if len(tmp) > 0:
                res[term] = set([mapping[x] for x in tmp])
        return res

def get_weights(label, eps):
    classes = np.array(np.array(label) * 10, dtype=np.int64)
    (unique, counts) = np.unique(classes, return_counts=True)
    no_counts = [x for x in range(unique[-1]) if x not in unique]
    for x in no_counts:
        counts = np.insert(counts, x, 0)
    weights = counts / counts.sum()
    weights = 1.0 / weights + eps
    weights = np.nan_to_num(weights, nan=0, neginf=0, posinf=0)
    weights = weights / weights.sum()

    return weights[classes]