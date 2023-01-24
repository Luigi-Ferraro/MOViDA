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

     

---------

### Environment set up
* Python libraries
     * networkx
     * numpy
     * pandas
     * torchmetrics
     * torchvision
