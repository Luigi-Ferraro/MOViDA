import torch
import numpy as np

def weighted_mse_loss(output, target, weights):
	classes = (target * 10).long()
	classes_output = (output * 10).long()
	classes_output = torch.clamp(classes_output, 0, weights.shape[0] - 1)
	return torch.sum(torch.max(weights[classes], weights[classes_output]) * (output - target) ** 2) / output.shape[0]

def mse_loss(output, target):
	return torch.sum((output - target) ** 2) / output.shape[0]

def mean_mse_loss(output, target):
	classes = (target * 10).long()
	idx_unique = classes.unique(sorted=True)
	errors = (output - target) ** 2
	return torch.stack([errors[classes == idx].mean() for idx in idx_unique]).mean()
 