import os
from pathlib import Path
import shutil

def create_dir(directory, rm_if_exists = False):
    dirpath = Path(directory)

    if not dirpath.exists():
        dirpath.mkdir(parents=True, exist_ok=True)
    elif dirpath.is_dir():
        if rm_if_exists:
            shutil.rmtree(dirpath)
            dirpath.mkdir(parents=True, exist_ok=True)
        

def create_train_info(filepath):
    with open(filepath, "w") as f:
        f.write("epoch,train_loss,train_metric,val_loss,val_metric,elapsed_time\n")	
        

def update_train_info(filepath, epoch,\
                        tot_loss_train, tot_metric_train, \
                        tot_loss_val, tot_metric_val, \
                        elapsed_time):
    with open(filepath, "a") as f:
        f.write(str(epoch) + "," + \
                str(tot_loss_train) + "," + str(tot_metric_train) + "," +\
                str(tot_loss_val) + "," + str(tot_metric_val) + "," +\
                str(elapsed_time))	
        

def create_test_info(filepath):
    path = Path(filepath)
    if not path.exists():
        with open(filepath, "w") as f:
            f.write("test_loss,test_metric,elapsed_time\n")	
        

def update_test_info(filepath,\
                        tot_loss_train, tot_metric_train, \
                        elapsed_time):
    with open(filepath, "a") as f:
        f.write(str(tot_loss_train) + "," + str(tot_metric_train) + "," +\
                str(elapsed_time))	