import torch
import torch.nn as nn
import torch.utils.data as du
from torch.utils.data.sampler import WeightedRandomSampler
import numpy as np

import time
import sys

from prepare_directories import *
from param_data import get_weights

class Model(nn.Module):
    def __init__(self, param_train, data_mv):
        
        super(Model, self).__init__()

        self.cuda           = param_train.cuda

        self.root           = data_mv.root
        self.ccl_hiddens    = param_train.ccl_hiddens
        self.drug_hiddens   = param_train.drug_hiddens
        self.final_hiddens  = param_train.final_hiddens

        self.mo_names       = data_mv.names_mo

        self.mo_gene2term   = data_mv.mo_gene2term

        self.genes_dim       = data_mv.genes_dim
        self.drugs_dim       = data_mv.drugs_dim

        self.construct_model(data_mv.hierarchy_graph)
        


    def construct_model(self, dG):
        self.contruct_VNN_input_layer()
        self.construct_VNN_graph(dG)

        self.construct_NN_drug()

        self.construct_NN_final()


    def contruct_VNN_input_layer(self):		
        for idx in range(len(self.mo_names)):
            for term, gene_set in self.mo_gene2term[idx].items():
                self.add_module(term + '_' + self.mo_names[idx] + '_gene_layer', nn.Linear(self.genes_dim[idx], len(gene_set)))


    def construct_VNN_graph(self, dG):

        self.VNN_layers = []   
        self.term_children_map = {}

        for term in dG.nodes():
            self.term_children_map[term] = []
            for child in dG.neighbors(term):
                self.term_children_map[term].append(child)

        while True:
            leaves = [n for n in dG.nodes() if dG.out_degree(n) == 0]

            if len(leaves) == 0:    
                break

            self.VNN_layers.append(leaves)

            for term in leaves:
				# input size will be #chilren + #genes directly annotated by the term
                input_size = 0
				# input size with one Pathway activity
                input_size += 1

                for child in self.term_children_map[term]:
                    input_size += self.ccl_hiddens
                    # input size with one Pathway activity for each child
                    input_size += 1

                for mo in self.mo_gene2term:
                    if term in mo:
                        input_size += len(mo[term])

                self.add_module(term+'_linear_layer', nn.Linear(input_size, self.ccl_hiddens))
                self.add_module(term+'_batchnorm_layer', nn.BatchNorm1d(self.ccl_hiddens))
            
            dG.remove_nodes_from(leaves)


    def construct_NN_drug(self):
        input_size = sum(self.drugs_dim)

        for i in range(len(self.drug_hiddens)):
            self.add_module('drug_linear_layer_' + str(i+1), nn.Linear(input_size, self.drug_hiddens[i]))
            self.add_module('drug_batchnorm_layer_' + str(i+1), nn.BatchNorm1d(self.drug_hiddens[i]))
            input_size = self.drug_hiddens[i]


    def construct_NN_final(self):
        final_input_size = self.ccl_hiddens + self.drug_hiddens[-1]
        self.add_module('final_linear_layer_1', nn.Linear(final_input_size, self.final_hiddens))
        self.add_module('final_batchnorm_layer_1', nn.BatchNorm1d(self.final_hiddens))
        self.add_module('final_linear_layer_2', nn.Linear(self.final_hiddens,1))


    def leakyReLU(self, x, slope = 0.1):
        return torch.relu(x) + slope * torch.min(torch.zeros(x.shape).cuda(self.cuda), x)

 
    def create_term_mask(self, term2gene_map, gene_dim):
        '''
        Build mask: matrix (nrows = number of relevant gene set, ncols = number all genes)
        elements of matrix are 1 if the corresponding gene is one of the relevant genes
        '''
        term_mask_map = {}
        for term, gene_set in term2gene_map.items():
            mask = torch.zeros(len(gene_set), gene_dim)
            for i, gene_id in enumerate(gene_set):
                mask[i, gene_id] = 1
            mask_gpu = torch.autograd.Variable(mask.cuda(self.cuda))
            term_mask_map[term] = mask_gpu
        return term_mask_map


    def initialize_model(self, learning_rate):

        self.to(self.cuda)

        self.optimizer = torch.optim.Adam(self.parameters(), lr=learning_rate, betas=(0.9, 0.99), eps=1e-05)
        self.optimizer.zero_grad()

        self.mo_term_masks = {moname : self.create_term_mask(gene2term, genedim) \
                        for moname, gene2term, genedim in zip(self.mo_names, self.mo_gene2term, self.genes_dim)}

        for name, param in self.named_parameters():
            term_name, mo_name = name.split('_')[:2]
            if '_gene_layer.weight' in name:
                param.data = torch.mul(param.data, self.mo_term_masks[mo_name][term_name]) * 0.1
            else:
                param.data = param.data * 0.1


    def reset_grad_mask(self):
        for name, param in self.named_parameters():
            if '_gene_layer.weight.weight' in name:
                term_name, mo_name = name.split('_')[:2]
                param.grad.data = torch.mul(param.grad.data, self.mo_term_masks[mo_name][term_name])


    def initialize_dataloader(self, dataset, eps_weights, batch_size, sampler_bool):
        sampler = None
        feature, label = dataset
        self.weights = get_weights(label, eps_weights)
        if sampler_bool:
            sampler = WeightedRandomSampler(self.weights, \
                                            label.shape[0])

        return du.DataLoader(du.TensorDataset(feature, label), batch_size=batch_size, \
                                                    shuffle=False, sampler = sampler)


    def extract_input(self, input_data, cell_features, drug_features, path_act, term2id_mapping):
        
        ccl_inputs      = input_data[:, 0].int()
        drug_inputs     = input_data[:, 1].int()

        ccl_sub         = [torch.autograd.Variable(torch.from_numpy(x[ccl_inputs,:]).float().cuda(self.cuda)) for x in cell_features]
        drug_sub        = [torch.from_numpy(x[drug_inputs,:]).float() for x in drug_features]
        drug_var        = torch.autograd.Variable(torch.cat(drug_sub,1)).cuda(self.cuda)
        path_act_var    = torch.from_numpy(path_act[ccl_inputs,:]).float().cuda(self.cuda)

        features        = (ccl_sub, drug_var, path_act_var, term2id_mapping)

        del ccl_sub, drug_var, path_act_var, term2id_mapping
        return features

    
    def model_train(self, train_data, val_data, ccl_features, drug_features, path_act, term2id_mapping, param_train):
        
        train_info_file = param_train.currexpdir + "train_info.csv"
        create_train_info(train_info_file)
        
        self.initialize_model(param_train.lr)
        train_loader    = self.initialize_dataloader(train_data, param_train.eps_weights,\
                                                     param_train.batchsize, param_train.sampler_bool)
        val_loader      = self.initialize_dataloader(val_data, param_train.eps_weights,\
                                                     param_train.batchsize, False)

        min_loss = None

        for epoch in range(param_train.epoch):
            epoch_start_time = time.time()

            self.train()
            tot_loss_train = []
            tot_metric_train = []

            for i, (inputdata, labels) in enumerate(train_loader):
                self.optimizer.zero_grad()  
                loss, metric, _ = self._step(inputdata, ccl_features, drug_features, path_act, term2id_mapping, labels, param_train)
                loss.backward()
                self.optimizer.step()

                tot_loss_train      += [loss.item()]
                tot_metric_train    += [metric.item()]

                self.reset_grad_mask()
            
            tot_loss_train      = np.mean(tot_loss_train)
            tot_metric_train    = np.mean(tot_metric_train)
            


            self.eval()
            tot_loss_val    = []
            tot_metric_val  = []
            with torch.no_grad():
                for i, (inputdata, labels) in enumerate(val_loader):
                    loss, metric, _ = self._step(inputdata, ccl_features, drug_features, path_act, term2id_mapping, labels, param_train)

                    tot_loss_val    += [loss.item()]
                    tot_metric_val  += [metric.item()]
                
                tot_loss_val    = np.mean(tot_loss_val)
                tot_metric_val  = np.mean(tot_metric_val)

            


            if min_loss is None or tot_loss_val < min_loss:
                min_loss = tot_loss_val
                torch.save(self, param_train.currexpdir + '/model_' + 'best' + '.pt')


            epoch_end_time = time.time()
            update_train_info(train_info_file, True, epoch,\
                                tot_loss_train, tot_metric_train, \
                                tot_loss_val, tot_metric_val, \
                                epoch_end_time - epoch_start_time)


    def model_test(self, test_data, ccl_features, drug_features, path_act, term2id_mapping, param_train):
        test_info_file = param_train.currexpdir + "test_info_file.csv"
        create_test_info(test_info_file)

        test_loader  = self.initialize_dataloader(test_data, param_train.eps_weights,\
                                                     param_train.batchsize, False)

        pred_tot         = None

        tot_loss_test    = []
        tot_metric_test  = []
        epoch_start_time = time.time()
        self.eval()
        with torch.no_grad():
            for i, (inputdata, labels) in enumerate(test_loader):
                loss, metric, pred = self._step(inputdata, ccl_features, drug_features, path_act, term2id_mapping, labels, param_train)
                
                tot_loss_test    += [loss.item()]
                tot_metric_test  += [metric.item()]
                
                tot_loss_test    = np.mean(tot_loss_test)
                tot_metric_test  = np.mean(tot_metric_test)

                if pred_tot is None:
                    pred_tot = pred.cpu().numpy()
                else:
                    pred_tot = np.concatenate([pred_tot, pred.cpu().numpy()])
            
            epoch_end_time = time.time()
            update_test_info(test_info_file, True, \
                                tot_loss_test, tot_metric_test, \
                                epoch_end_time - epoch_start_time)

        return pred_tot


    def _step(self, inputdata, ccl_features, drug_features, path_act, term2id_mapping, labels, param_train):
        features = self.extract_input(inputdata, ccl_features, drug_features, path_act, term2id_mapping)
            

        pred = self.forward(features)
        del features

        if labels is not None:
            cuda_labels = torch.autograd.Variable(labels.cuda(self.cuda))
            loss = param_train.loss_func(pred, cuda_labels, torch.from_numpy(self.weights).cuda(self.cuda))
            metric = param_train.metric_func(pred, cuda_labels)

            del cuda_labels

            return loss, metric, pred
        else:
            return pred


    def forward_gene_layer(self, mo_ccl_input):
        self.term_gene_out_mo = {mo_name : {} for mo_name in self.mo_names}
        for idx_mo_name,mo in zip(range(len(self.mo_names)), self.mo_gene2term):
            for term in mo.keys():
                self.term_gene_out_mo[self.mo_names[idx_mo_name]][term] = \
                    self._modules[term + "_" + self.mo_names[idx_mo_name] + '_gene_layer'](mo_ccl_input[idx_mo_name])


    def forward_VNN_layers(self, path_act_input, term2id_mapping):
        self.term_VNN_out = {}

        for i, layer in enumerate(self.VNN_layers):

            for term in layer:
                children_term = self.term_children_map[term]
                ids_terms = [term2id_mapping[t] for t in [term] + children_term]
                path_act_term = path_act_input[:,ids_terms].reshape(-1, len(children_term) + 1)
                path_act_term = torch.autograd.Variable(path_act_term.type(torch.FloatTensor).cuda(self.cuda))

                child_input_list = []

                for child in children_term:
                    child_input_list.append(self.term_VNN_out[child])
                
                for mo in self.mo_names:
                    if term in self.term_gene_out_mo[mo]:
                        child_input_list.append(self.term_gene_out_mo[mo][term])

                child_input_list.append(path_act_term)
                del path_act_term

                child_input_list = torch.cat(child_input_list,1)

                term_NN_out = self._modules[term+'_linear_layer'](child_input_list)		
                del child_input_list		

                term_NN_out = torch.tanh(term_NN_out)
                term_NN_out = self._modules[term+'_batchnorm_layer'](term_NN_out)
                self.term_VNN_out[term] = term_NN_out
        del self.term_VNN_out
        return term_NN_out

        
    def forward_drugs_layers(self, drugs_input):
        drug_out = drugs_input

        #self.term_NN_drug_out = {}
        for i in range(1, len(self.drug_hiddens)+1, 1):
            drug_out = self._modules['drug_batchnorm_layer_'+str(i)]( torch.tanh(self._modules['drug_linear_layer_' + str(i)](drug_out)))
            #self.term_NN_drug_out['drug_'+str(i)] = drug_out
        
        return drug_out


    def forward_final_layers(self, ccl_out, drug_out):
        self.term_NN_final_out = {}
        final_input = torch.cat((ccl_out, drug_out), 1)

        out = self._modules['final_batchnorm_layer_1'](torch.tanh(self._modules['final_linear_layer_1'](final_input)))
        self.term_NN_final_out["final_1"] = out

        out = self.leakyReLU(self._modules['final_linear_layer_2'](out))
        self.term_NN_final_out['final'] = out
        return out


    def forward(self, x):
        mo_ccl_input, drugs_input, path_act_input, term2id_mapping = x

        self.forward_gene_layer(mo_ccl_input)
        del mo_ccl_input
        
        ccl_out = self.forward_VNN_layers(path_act_input, term2id_mapping)
        del path_act_input, term2id_mapping

        drug_out = self.forward_drugs_layers(drugs_input)
        del drugs_input

        return self.forward_final_layers(ccl_out, drug_out)
