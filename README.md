## MOViDA

## Introduction
Our tool makes drug sensitivity prediction for cancer cell lines using a XAI model based on gene ontology hierarchy. We use multi-omics data, including mutations, amplifications, deletions and pathway activity to make predictions. This cutting-edge approach allows for a more accurate prediction of drug sensitivity by taking into account the specific genetic makeup of the cancer cell line. By utilizing gene ontology hierarchy and multi-omics data, our tool provides a powerful resource for researchers and clinicians looking to improve treatment outcomes for patients with cancer.


---------

### Environment set up

* Hardware: GPU with CUDA >= 11
* Software:
    * Python >= 3.7
    * [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
    * PyTorch >= 11.3 (verify on [Pytorch](https://pytorch.org/get-started/locally/) the corrected version for your pc)
        ```angular2
        conda install pytorch torchvision torchaudio pytorch-cuda=11.7 -c pytorch -c nvidia
        ```
    * Python libraries
         * networkx
         * numpy
         * pandas
         * torchmetrics
         * torchvision



---------

### Files

* code
   * main.py : starts training or testing
   * param_data.py : contains a class that represents the data in input 
   * param_train.py : contains a class that represents the training parameters
   * prepare_directories.py : creates directory and file for the experiment
   * models : directory that contains classes which encode all the parameters and structure of the models, as well as functions needed for training them
      * MOViDA.py : class for drug sensitivity prediction
      * MOViDA_synergy.py : class for synergistic drug combination prediction
   * parameters : directory that contains .txt files for data and parameters information, see below


---------

### Parameters
A parameters file is a .txt file that contains a Python dictionary with all the necessary information for data and training. An example of such a file would be:
   ```
   {
       "input_dir"         : "data", 
       "names_mo"          : ["mut", "amp", "del"], 
       ...
   }
   ```

* data parameters : This is a python dictionary that contains the input data needed to launch a training process. The dictionary contains the following key-value pairs :
   * "input_dir": This is the path to the directory where the input data files are located.
   * "names_mo": This is a list of strings that contains the names of the different types of genetic modifications (e.g. "mut", "amp", "del")
   * "drug2id": This is the path to the file that contains the mapping from drug names to drug IDs.
   * "cell2id": This is the path to the file that contains the mapping from cell line names to cell line IDs.
   * "mo_gene2id": This is a list of strings that contains the paths to the files that map the gene names to the genetic modification IDs for each type of modification.
   * "train_file": This is the path to the file that contains the training set of drug-cell line pairs.
   * "val_file": This is the path to the file that contains the validation set of drug-cell line pairs.
   * "test_file": This is the path to the file that contains the test set of drug-cell line pairs.
   * "all_file": This is the path to the file that contains all the drug-cell line pairs.
   * "cell2features": This is a list of strings that contains the paths to the files that contain the features of the cell lines for each type of modification.
   * "names_drfeat": This is a list of strings that contains the names of the different types of drug features (e.g. "pc", "vs").
   * "drug2features": This is a list of strings that contains the paths to the files that contain the features of the drugs for each type of feature.
   * "pathway_act_file": This is the path to the file that contains the pathway activity data.
   * "ontology": This is the path to the file that contains the ontology data.
   * "mo_gene2ontology": This is a list of strings that contains the paths to the files that map the genes to the ontology for each type of modification.
   * "synergy_bool": This is a boolean value that indicates whether the training process should include synergy information (i.e. whether the output should predict synergy scores or not).

* training parameters : This is a python dictionary that contains the training parameters needed to launch a training process :
   * "epoch": This is the number of training iterations.
   * "lr": This is the learning rate used for training.
   * "batchsize": This is the number of samples in a batch used for training.
   * "loss_func": This is the name of the loss function used for training.
   * "metric_func": This is the name of the metric used to evaluate the performance of the model.
   * "classifier": This is a boolean value that indicates whether the training process is for classification or regression.
   * "wloss_bool": This is a boolean value that indicates whether the loss function is weighted or not.
   * "focal_bool": This is a boolean value that indicates whether the Focal Loss is used or not.
   * "f_alpha": This is the value used for Focal Loss alpha.
   * "f_gamma": This is the value used for Focal Loss gamma.
   * "sampler_bool": This is a boolean value that indicates whether the data should be oversampled or not.
   * "ccl_hiddens": This is the number of nodes in hidden layers in the cell line encoder.
   * "drug_hiddens": This is a string that contains the number of hidden layers in the drug encoder and how many nodes they have.
   * "final_hiddens": This is the number of nodes in hidden layers in the final network.
   * "max_val": This is the maximum value used for the target variable.
   * "eps_weights": This is the parameter epsilon using to balance the weights.
   * "expdir": This is the directory where all the experiments results will be stored.
   * "currexpdir": This is the directory where the current experiment results will be stored, it will be in "expdir".
   * "result_file": This is the path to the file where the test set predictions will be stored.
   * "result_all_file": This is the path to the file where the all set predictions will be stored.
   * "model_name": This is the name of python file storing the model used (without the file extension .py)
   * "load_model": This is the path to the pre-trained model if any.




---------

### Start training
* prepare dictionary files for data and parameter
* execute this command whit the correct data_dict and param_dict
   ```
   python code/main.py -data_dict data_dict -param_dict param_dict -cuda 0
   ```

