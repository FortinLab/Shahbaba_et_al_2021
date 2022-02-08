from typing import List, Any
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier
import numpy as np
import time
import os

import torch
import torch.nn as nn
import torch.optim as optim
from torch.optim.optimizer import Optimizer


class FullyConnectedModel(nn.Module):
    def _forward_unimplemented(self, *input_val: Any) -> None:
        pass

    def __init__(self, input_dim: int, layer_dims: List[int], layer_batch_norms: List[bool], layer_acts: List[str]):
        super().__init__()
        self.layers: nn.ModuleList[nn.ModuleList] = nn.ModuleList()

        # layers
        for layer_dim, batch_norm, act in zip(layer_dims, layer_batch_norms, layer_acts):
            module_list = nn.ModuleList()

            # linear
            module_list.append(nn.Linear(input_dim, layer_dim))

            # batch norm
            if batch_norm:
                module_list.append(nn.BatchNorm1d(layer_dim))

            # activation
            act = act.upper()
            if act == "RELU":
                module_list.append(nn.ReLU())
            elif act == "SIGMOID":
                module_list.append(nn.Sigmoid())
            elif act != "LINEAR":
                raise ValueError("Un-defined activation type %s" % act)

            self.layers.append(module_list)

            input_dim = layer_dim

    def forward(self, x):
        x = x.float()

        module_list: nn.ModuleList
        for module_list in self.layers:
            for module in module_list:
                x = module(x)

        return x


def get_enc_model(input_dim: int, enc_dim: int) -> nn.Module:
    h_dim = 500
    n_layers = 2

    nnet_enc = FullyConnectedModel(input_dim, [h_dim] * n_layers + [enc_dim], [False] * n_layers + [False],
                                   ["RELU"] * n_layers + ["LINEAR"])

    return nnet_enc


def get_dec_model(input_dim: int, enc_dim: int) -> nn.Module:
    h_dim = 500
    n_layers = 2

    nnet_dec = FullyConnectedModel(enc_dim, [h_dim] * n_layers + [input_dim], [False] * n_layers + [False],
                                   ["RELU"] * n_layers + ["LINEAR"])

    return nnet_dec


def train_ae(data: np.ndarray, enc_dim: int, model_save_loc: str, device: torch.device):
    if not os.path.exists(model_save_loc):
        os.makedirs(model_save_loc)

    input_dim: int = int(data.shape[1])

    # define training hyperparams
    lr_init = 0.1
    lr_decay = 0.9999993
    mom_init = 0.5
    mom_final = 0.9
    mom_change_steps = 10000
    display_itrs = 1000
    # num_itrs = 200
    num_itrs = 20000
    batch_size = 1000

    # get neural network
    nnet_enc = get_enc_model(input_dim, enc_dim).to(device)
    nnet_dec = get_dec_model(input_dim, enc_dim).to(device)

    nnet_enc.train()
    nnet_dec.train()

    criterion = nn.MSELoss()
    optimizer: Optimizer = optim.SGD(list(nnet_enc.parameters()) + list(nnet_dec.parameters()), lr=lr_init,
                                     momentum=mom_init)
    # optimizer: Optimizer = optim.Adam(list(nnet_enc.parameters()) + list(nnet_dec.parameters()), lr=0.001)

    display_start_time = time.time()
    for train_itr in range(num_itrs):
        # zero the parameter gradients
        optimizer.zero_grad()
        lr_itr: float = lr_init * (lr_decay ** train_itr)
        mom_itr = mom_init + (mom_final - mom_init) * min(train_itr / mom_change_steps, 1.0)

        for param_group in optimizer.param_groups:
            param_group['momentum'] = mom_itr
            param_group['lr'] = lr_itr

        # get batch data
        rand_idxs = np.random.choice(data.shape[0], batch_size, replace=True)
        data_batch = torch.tensor(data[rand_idxs], device=device).float()

        # get autoencoding
        data_batch_enc = nnet_enc(data_batch)
        data_batch_dec = nnet_dec(data_batch_enc)

        # get loss
        loss = criterion(data_batch, data_batch_dec)

        # backwards
        loss.backward()

        # optimizer step
        optimizer.step()

        if train_itr % display_itrs == 0:
            display_elapsed_time = time.time() - display_start_time

            print("Itr: %i, lr: %f, mom: %.2f, cost: %s, time: %.2f" % (
                train_itr, lr_itr, mom_itr, loss.item(), display_elapsed_time))
            display_start_time = time.time()

    print("Saving network to %s" % model_save_loc)
    torch.save(nnet_enc.state_dict(), "%s/model_state_dict_enc.pt" % model_save_loc)
    torch.save(nnet_dec.state_dict(), "%s/model_state_dict_dec.pt" % model_save_loc)

    return nnet_enc


# pytorch device
def get_device() -> torch.device:
    if ('CUDA_VISIBLE_DEVICES' in os.environ) and torch.cuda.is_available():
        device = torch.device("cuda:%i" % 0)
    else:
        torch.set_num_threads(1)
        device: torch.device = torch.device("cpu")

    return device


def pca_enc_fn(data, dim):
    pca = PCA(n_components=dim)
    pca.fit(data)

    def enc_fn(x):
        x = x.reshape((x.shape[0], np.prod(x.shape[1:])))
        enc_r = pca.transform(x)
        return enc_r

    return enc_fn


# Classification methods
def knn_fn(data_x, data_y, knn):
    classifier = KNeighborsClassifier(n_neighbors=knn)
    classifier.fit(data_x, data_y)

    return classifier.predict
