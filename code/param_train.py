class ParamTrainMOViDA:
    epoch               = None
    lr                  = None
    batchsize           = None
    loss_func           = None
    metric_func         = None
    sampler_bool        = None
    expdir              = None
    result_file         = None
    result_all_file     = None
    cuda                = None
    ccl_hiddens         = None
    drug_hiddens        = None
    final_hiddens       = None
    max_val             = None
    eps_weights         = None
    model_name          = None

    def __init__(self, cuda, 
                epoch, lr, batchsize, 
                loss_func, metric_func, sampler_bool, 
                ccl_hiddens, drug_hiddens, final_hiddens, 
                max_val, eps_weights,
                expdir, currexpdir, result_file, result_all_file,
                model_name, load_model):

        self.cuda               = cuda

        self.epoch              = epoch
        self.lr                 = lr
        self.batchsize          = batchsize

        self.loss_func          = loss_func
        self.metric_func        = metric_func

        self.sampler_bool       = sampler_bool

        self.ccl_hiddens        = ccl_hiddens
        self.drug_hiddens       = list(map(int, drug_hiddens.split(',')))
        self.final_hiddens      = final_hiddens

        self.max_val            = max_val
        self.eps_weights        = eps_weights

        self.expdir             = expdir + "/"
        self.currexpdir         = self.expdir + currexpdir + "/"
        self.result_file        = self.currexpdir + result_file
        self.result_all_file    = self.currexpdir + result_all_file

        self.model_name         = model_name

        self.load_model         = load_model