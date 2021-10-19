import numpy as np
import argparse

from utils import data_utils
from utils import ml_utils
import pickle


def shuffle_labels(labels):
    labelName_first = list(labels.keys())[0]
    numLabels = labels[labelName_first].shape[0]

    randIdxs = np.random.choice(numLabels, numLabels, replace=False)

    for key in labels.keys():
        labels[key] = labels[key][randIdxs]

    return labels


def classification(dataset_name, time_window, inp_len, num_shuffles):
    knn = 2

    data = data_utils.get_dataset(dataset_name, time_window, inp_len)

    accs = []

    encs = pickle.load(open("data/%s/encoding_%s_%s" % (dataset_name, time_window, inp_len), "rb"))

    """
    labelsSelects = [
                    dict(),
                    dict(),
                    {"odor":[1,2,3],"isCorrect":[1]},
                    {"isInSeq":[1],"odor":[1,2,3],"isCorrect":[1]},
                    {"odor":[1,2,3],"isCorrect":[0]},
                    {"isInSeq":[1],"odor":[1,2,3],"isCorrect":[0]},

                    dict()
                   ]

    accMeasures = ["isInSeq","odor","isInSeq","odor","isInSeq","odor","isCorrect"]
    """
    labelsSelects = [
        {"odor": [1, 2, 3], "isCorrect": [1]},
        {"isInSeq": [1], "odor": [1, 2, 3], "isCorrect": [1]},
        dict()
    ]

    accMeasures = ["isInSeq", "odor", "isCorrect"]

    for expNum in range(num_shuffles + 1):
        print("Run: (%i/%i)" % (expNum, num_shuffles))

        acc = []
        for encs_t, data_t in zip(encs, data):
            acc_t = dict()
            acc_t["startTime"] = min(encs_t["times"])

            labelsTrain_t = data_t["labelsTrain"]
            labelsTest_t = data_t["labelsTest"]
            if expNum > 0:
                # Shuffle data
                labelsTrain_t = shuffle_labels(labelsTrain_t)
                labelsTest_t = shuffle_labels(labelsTest_t)

            # Make labels
            for labelsSelect, accMeasure in zip(labelsSelects, accMeasures):
                # Train model
                _, input_select_train, labels_select_train, name = data_utils.get_selected_data(labelsSelect,
                                                                                                accMeasure,
                                                                                                encs_t["inputTrain"],
                                                                                                labelsTrain_t)
                if input_select_train.shape[0] == 0:
                    print("No training examples for %s for accMeasure %s on trial %i, startTime %.3f, "
                          "assigning NAN" % (dataset_name, accMeasure, expNum, acc_t["startTime"]))
                    acc_t[name] = float('nan')
                    continue
                predModel = ml_utils.knn_fn(input_select_train, labels_select_train, knn)

                # Test model
                _, input_select_test, labels_select_test, name = data_utils.get_selected_data(labelsSelect, accMeasure,
                                                                                              encs_t["inputTest"],
                                                                                              labelsTest_t)
                if input_select_train.shape[0] == 0:
                    print("No testing examples for %s for accMeasure %s on trial %i, startTime %.3f, "
                          "assigning NAN" % (dataset_name, accMeasure, expNum, acc_t["startTime"]))
                    acc_t[name] = float('nan')
                    continue

                labels_test_pred_t = predModel(input_select_test)
                acc_test = 100.0 * np.sum(labels_select_test == labels_test_pred_t) / float(input_select_test.shape[0])

                acc_t[name] = acc_test

            acc.append(acc_t)

        accs.append(acc)

    pickle.dump(accs, open("data/%s/accs_%s_%s" % (dataset_name, time_window, inp_len), "wb"))

    return accs


def main():
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--timewindow', type=float, required=True, help="Time window of plots")
    parser.add_argument('--inpLen', type=int, required=True, help="Input length")
    parser.add_argument('--num_shuffles', type=int, default=0, help="Number of random shuffles")
    parser.add_argument('--datasets', type=str, required=True, help="Names of datasets separated by a comma")

    args = parser.parse_args()

    dataset_names = args.datasets.split(",")
    for dataset_name in dataset_names:
        print("DATASET: %s" % dataset_name)
        classification(dataset_name, args.timewindow, args.inpLen, args.num_shuffles)


if __name__ == "__main__":
    main()
