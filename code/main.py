#%%
import os
os.chdir('/storage/qnap_home/ferraro/Progetti/MOViDA')
os.getcwd()

#%%
import numpy as np
import torch
import importlib

# %%
from param_data             import *
from param_train            import *
from prepare_directories    import *
from losses                 import *


#%%
# Start reproducibility
import random
torch.manual_seed(0)
np.random.seed(0)
random.seed(0)

torch.cuda.manual_seed(0)
torch.backends.cudnn.deterministic = True  # Note that this Deterministic mode can have a performance impact
torch.backends.cudnn.benchmark = False
# End reproducibility


# %%
param_train = ParamTrainMOViDA(
                cuda                = 0, 
                epoch               = 2, 
                lr                  = 0.001, 
                batchsize           = 10000, 
                loss_func           = weighted_mse_loss,
                metric_func         = mean_mse_loss,
                sampler_bool        = True,
                ccl_hiddens         = 6,
                drug_hiddens        = "100,50,6", 
                final_hiddens       = 6, 
                max_val             = 1.2,
                eps_weights         = 80,
                expdir              = "experiments", 
                currexpdir          = "prova",
                result_file         = "pred_test.txt",
                result_all_file     = "pred_all.txt",
                model_name          = "MOViDA",
                load_model          = None
)


# %%
data_mv = DataMOViDA(
                input_dir           = "data_tmp", 
                drug2id             = "drug2ind.txt", 
                cell2id             = "cell2ind.txt", 
                mo_gene2id          = ["gene2ind2mut.txt", "gene2ind2amp.txt", "gene2ind2del.txt"],
                train_file          = "drugcell_train.txt", 
                val_file            = "drugcell_val.txt", 
                test_file           = "drugcell_test.txt",
                all_file            = "drugcell_all.txt",
                cell2features       = ["cell2mutation.txt", "cell2amp.txt", "cell2del.txt"], 
                drug2features       = ["drug2pc_fingerprint.txt", "drugfeatnew_all.txt"], 
                pathway_act_file    = "GO_enrich.tsv",
                ontology            = "ont_tree.txt", 
                mo_gene2ontology    = ["ont_mut.txt", "ont_amp.txt", "ont_del.txt"]
)

def train():
    #%%
    module = importlib.import_module("models." + param_train.model_name)
    # %%
    create_dir(param_train.expdir)
    create_dir(param_train.currexpdir, rm_if_exists = True)

    print("Create Model")

    # %%
    model = module.Model(param_train, data_mv)

    print("Start training")

    # %%
    model.model_train(data_mv.train_data, data_mv.val_data, data_mv.ccl_features, data_mv.drug_features, param_train)

    print("Start testing")

    # %%
    pred = model.model_test(data_mv.test_data, data_mv.ccl_features, data_mv.drug_features, param_train)
    np.savetxt(param_train.result_file, pred,'%.4e')
    
    pred = model.model_test(data_mv.all_data, data_mv.ccl_features, data_mv.drug_features, param_train)
    np.savetxt(param_train.result_all_file, pred,'%.4e')

    print("END")

    return pred

def test():
    
    print("Start testing")

    model = torch.load(param_train.load_model, map_location='cuda:%d' % param_train.cuda)
    pred = model.model_test(data_mv.test_data, data_mv.ccl_features, data_mv.drug_features, param_train)
    np.savetxt(param_train.result_file, pred,'%.4e')

    print("END")
    
    return pred


if __name__ == "__main__":
    #args = [dict_data, dict_param]

    if param_train.load_model is None:
        train()
    else:
        test()
