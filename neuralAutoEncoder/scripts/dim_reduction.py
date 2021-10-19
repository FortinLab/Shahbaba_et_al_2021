from utils import data_utils
from utils import ml_utils

import torch
import torch.nn as nn

import numpy as np
import argparse

import pickle


def dim_reduction(dataset_name, time_window, inp_len):
    enc_dim = 2

    data = data_utils.get_dataset(dataset_name, time_window, inp_len)

    encs = []

    method: str = "NNET"

    device: torch.device = ml_utils.get_device()
    if method.upper() == "NNET":
        print(device)

    for data_t in data:
        start_time = min(data_t["times"])
        end_time = max(data_t["times"])
        print("Start time %s, End Time %s" % (start_time, end_time))

        encs_t = dict()
        encs_t["times"] = data_t["times"]

        # Get input data
        input_train_t = data_t["inputTrain"]
        input_test_t = data_t["inputTest"]
        train_test_t = np.concatenate([input_train_t, input_test_t], axis=0)
        print("%i training examples" % train_test_t.shape[0])

        # Train enc_dim reduction method for time window
        if method.upper() == "NNET":
            # train nnet
            model_save_loc = "savedModels/%s_start%.3f_end%.3f_inpLen%i" % (dataset_name, start_time, end_time, inp_len)
            nnet_enc: nn.Module = ml_utils.train_ae(train_test_t, enc_dim, model_save_loc, device)
            nnet_enc.eval()

            # get encoding
            input_train_t_tens = torch.tensor(input_train_t, device=device).float()
            input_test_t_tens = torch.tensor(input_test_t, device=device).float()

            encs_t["inputTrain"] = nnet_enc(input_train_t_tens).cpu().data.numpy()
            encs_t["inputTest"] = nnet_enc(input_test_t_tens).cpu().data.numpy()

        elif method.upper() == "PCA":
            enc_fn = ml_utils.pca_enc_fn(train_test_t, enc_dim)

            # Do dimensionality reduction for data in time window
            encs_t["inputTrain"] = enc_fn(input_train_t)
            encs_t["inputTest"] = enc_fn(input_test_t)
        else:
            raise ValueError("Unknown method %s" % method)

        encs.append(encs_t)

    pickle.dump(encs, open("data/%s/encoding_%s_%s" % (dataset_name, time_window, inp_len), "wb"), protocol=1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--timewindow', type=float, required=True, help="Time window of plots")
    parser.add_argument('--inpLen', type=int, required=True, help="Input length")
    parser.add_argument('--datasets', type=str, required=True, help="Names of datasets separated by a comma")

    args = parser.parse_args()

    datasetNames = args.datasets.split(",")
    for dataset_name in datasetNames:
        print("DATASET: %s" % dataset_name)
        dim_reduction(dataset_name, args.timewindow, args.inpLen)


if __name__ == "__main__":
    main()
