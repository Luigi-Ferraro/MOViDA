import torch
import torch.nn as nn
import torch.utils.data as du
from torch.utils.data.sampler import WeightedRandomSampler
import numpy as np

import time

from prepare_directories import *
from param_data import get_weights

class Model(nn.Module):
    def __init__(self, param_train, data_mv):

        self.cuda           = param_train.cuda
        self.learning_rate  = param_train.lr
 
    def create_term_mask(self, term_direct_gene_map, gene_dim):
        '''
        Build mask: matrix (nrows = number of relevant gene set, ncols = number all genes)
        elements of matrix are 1 if the corresponding gene is one of the relevant genes
        '''
        term_mask_map = {}
        for term, gene_set in term_direct_gene_map.items():
            mask = torch.zeros(len(gene_set), gene_dim)
            for i, gene_id in enumerate(gene_set):
                mask[i, gene_id] = 1
            mask_gpu = torch.autograd.Variable(mask.cuda(self.cuda))
            term_mask_map[term] = mask_gpu
        return term_mask_map


    def initialize_model(self, data_mv):

        self.to(self.cuda)

        self.optimizer = torch.optim.Adam(self.parameters(), lr=self.learning_rate, betas=(0.9, 0.99), eps=1e-05)
        self.optimizer.zero_grad()

        mo_term_masks = [self.create_term_mask(gene2term, genedim) for gene2term, genedim in zip(data_mv.mo_gene2term, data_mv.genes_dim)]

        for name, param in self.named_parameters():
            term_name = name.split('_')[0]
            if '_direct_gene_mul_layer.weight' in name:
                param.data = torch.mul(param.data, term_mask_map_mul[term_name]) * 0.1
            elif '_direct_gene_gexpr_layer.weight' in name:
                param.data = torch.mul(param.data, term_mask_map_gexpr[term_name]) * 0.1
            else:
                param.data = param.data * 0.1


    def reset_grad_mask(self):
        for name, param in self.named_parameters():
            if '_direct_gene_layer.weight' not in name:
                continue
            term_name = name.split('_')[0]
            #print name, param.grad.data.size(), term_mask_map[term_name].size()
            if '_direct_gene_mul_layer.weight' in name:
                param.grad.data = torch.mul(param.grad.data, self.term_mask_map_mul[term_name])
            elif '_direct_gene_gexpr_layer.weight' in name:
                param.grad.data = torch.mul(param.grad.data, self.term_mask_map_gexpr[term_name])


    def initialize_dataloader(self, dataset, eps_weights, batch_size, sampler_bool):
        sampler = None
        feature, label = dataset
        if sampler_bool:
            sampler = WeightedRandomSampler(get_weights(label, eps_weights).squeeze(1), \
                                            label.shape[0])

        return du.DataLoader(du.TensorDataset(feature, label), batch_size=batch_size, \
                                                    shuffle=False, sampler = sampler)


    def extract_input(self, input_data, cell_features, drug_features):
        
        ccl_inputs  = input_data[:, 0]
        drug_inputs = input_data[:, 0]

        ccl_sub     = [torch.from_numpy(x[ccl_inputs]).float() for x in cell_features]
        drug_sub    = [torch.from_numpy(x[drug_inputs]).float() for x in drug_features]

        features    = (ccl_sub, drug_sub)
        return features

    
    def model_train(self, train_data, val_data, ccl_features, drug_features, param_train):
        
        train_info_file = param_train.expdir + "train_info.csv"
        create_train_info(train_info_file)
        
        self.initialize_model()
        train_loader    = self.initialize_dataloader(train_data, param_train.eps_weights,\
                                                     param_train.batch_size, param_train.sampler_bool)
        val_loader      = self.initialize_dataloader(val_data, param_train.eps_weights,\
                                                     param_train.batch_size, False)

        min_loss = None

        for epoch in range(param_train.epoch):
            epoch_start_time = time.time()

            self.train()
            tot_loss_train = 0
            tot_metric_train = 0

            for i, (inputdata, labels) in enumerate(train_loader):
                self.optimizer.zero_grad()  
                loss, metric, _ = self._step(inputdata, ccl_features, drug_features, labels, param_train)
                loss.backward()
                self.optimizer.step()

                tot_loss_train      += loss.item()
                tot_metric_train    += metric

                self.reset_grad_mask()
            
            tot_loss_train      = torch.mean(tot_loss_train)
            tot_metric_train    = torch.mean(tot_metric_train)
            


            self.eval()
            tot_loss_val    = 0
            tot_metric_val  = 0

            for i, (inputdata, labels) in enumerate(val_loader):
                loss, metric, _ = self._step(inputdata, ccl_features, drug_features, labels, param_train)

                tot_loss_val    += loss.item()
                tot_metric_val  += metric
            
            tot_loss_val    = torch.mean(tot_loss_val)
            tot_metric_val  = torch.mean(tot_metric_val)

            


            if min_loss is None or tot_loss_val < min_loss:
                min_loss = tot_loss_val
                torch.save(self, param_train.currexpdir + '/model_' + 'best' + '.pt')


            epoch_end_time = time.time()
            update_train_info(train_info_file, epoch,\
                                tot_loss_train, tot_metric_train, \
                                tot_loss_val, tot_metric_val, \
                                epoch_end_time - epoch_start_time)


                

    def model_test(self, test_data, ccl_features, drug_features, param_train):
        test_info_file = param_train.expdir + "test_info_file.csv"
        create_test_info(test_info_file)

        test_loader  = self.initialize_dataloader(test_data, param_train.eps_weights,\
                                                     param_train.batch_size, False)

        pred_tot         = None

        tot_loss_test    = 0
        tot_metric_test  = 0
        epoch_start_time = time.time()
        self.eval()
        for i, (inputdata, labels) in enumerate(test_loader):
            loss, metric, pred = self._step(inputdata, ccl_features, drug_features, labels, param_train)
            
            tot_loss_test    += loss.item()
            tot_metric_test  += metric
            
            tot_loss_test    = torch.mean(tot_loss_test)
            tot_metric_test  = torch.mean(tot_metric_test)

            if pred_tot is not None:
                pred_tot = pred
            else:
                pred_tot = torch.cat(pred, pred_tot)
        
        epoch_end_time = time.time()
        update_train_info(test_info_file,\
                            tot_loss_test, tot_metric_test, \
                            epoch_end_time - epoch_start_time)

        return pred_tot.cpu().numpy()




    def _step(self, inputdata, ccl_features, drug_features, labels, param_train):
        features = self.extract_input(inputdata, ccl_features, drug_features)

        cuda_features = torch.autograd.Variable(features.cuda(self.cuda))
        if labels is not None:
            cuda_labels = torch.autograd.Variable(labels.cuda(self.cuda))

        pred = self.forward(cuda_features)

        if labels is not None:
            loss = param_train.loss_func(pred, cuda_labels, self.train_weights)
            metric = param_train.metric_func(pred, cuda_labels)

            return loss, metric, pred
        else:
            return pred
