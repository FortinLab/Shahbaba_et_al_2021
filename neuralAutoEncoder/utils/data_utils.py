import numpy as np


def flatten_data(input_data, labels, times, inp_len):
    dataFlat = []

    labelsFlat = dict()
    labelNames = labels.keys()
    for labelName in labelNames:
        labelsFlat[labelName] = []

    timesFlat = []
    encIdxs = range(0, input_data.shape[1] - inp_len + 1)
    for startIdx in encIdxs:
        endIdx = startIdx + inp_len
        data_slice = input_data[:, startIdx:endIdx, :]  # Get slice
        times_slice = times[startIdx:endIdx]  # Get slice

        dataFlat.append(data_slice.copy())
        for labelName in labelNames:
            labelsFlat[labelName].append(labels[labelName].copy())

        timesFlat.append(times_slice.copy())

    dataFlat = np.concatenate(dataFlat)
    for labelName in labelNames:
        labelsFlat[labelName] = np.concatenate(labelsFlat[labelName])

    timesFlat = np.concatenate(timesFlat)

    dataFlat = dataFlat.reshape((dataFlat.shape[0], dataFlat.shape[1] * dataFlat.shape[2]))

    return dataFlat, labelsFlat, timesFlat


def get_dataset(dataset_name, time_window, inp_len):
    dataset = load_raw_data(dataset_name)

    startTime = np.min(dataset["times"])
    maxTime = np.max(dataset["times"])

    dataset_windows = []

    while startTime < maxTime:
        endTime = startTime + time_window
        dataset_w = dict()

        idxs = (dataset["times"] >= startTime) & (dataset["times"] <= endTime)

        times_w = dataset["times"][idxs]
        inputTrain_w = dataset["inputTrain"][:, idxs, :]

        dataset_w["inputTrain"], dataset_w["labelsTrain"], dataset_w["times"] = flatten_data(inputTrain_w,
                                                                                             dataset["labelsTrain"],
                                                                                             times_w, inp_len)

        inputTest_w = dataset["inputTest"][:, idxs, :]
        dataset_w["inputTest"], dataset_w["labelsTest"], dataset_w["times"] = flatten_data(inputTest_w,
                                                                                           dataset["labelsTest"],
                                                                                           times_w, inp_len)

        dataset_windows.append(dataset_w)

        startTime = endTime

    return dataset_windows


def load_raw_data(rat_name):
    # Get dataset
    spikeTrain = np.load("data/%s/%s_PokeOut_spike_data_binned.npy" % (rat_name, rat_name))
    spikeTrain_labels = np.load("data/%s/%s_trial_info.npy" % (rat_name, rat_name))
    # spikeTrain = np.load("data_old/superchris/spike_data_binned.npy")
    # spikeTrain_labels = np.load("data_old/superchris/trial_info.npy")

    # Preprocess data
    times = np.arange(-2, 2, 0.01)
    spikeTrain = spikeTrain * 10
    spikeTrain = np.transpose(spikeTrain, [0, 2, 1])

    spikeTrain_labels[:, 2] = spikeTrain_labels[:, 2] - 1
    spikeTrain_labels[:, 3] = spikeTrain_labels[:, 3] - 1

    # print("Original: %i examples, %i timepoints, %i inputs per timepoint" % (spikeTrain.shape[0],
    # spikeTrain.shape[1],spikeTrain.shape[2]))

    # Remove instances where mouse was incorrect
    """  TODO, remove incorrect?
    isCorrect = spikeTrain_labels[:,0] == 1

    spikeTrain = spikeTrain[isCorrect]
    spikeTrain_labels = spikeTrain_labels[isCorrect]
    """

    # print("After: %i examples, %i timepoints, %i inputs per timepoint" %
    # (spikeTrain.shape[0],spikeTrain.shape[1],spikeTrain.shape[2]))

    dataset = dict()
    numTrials = spikeTrain.shape[0]
    trialID = np.arange(0, numTrials, 1)

    trainEndIdx = int(np.floor(numTrials * 0.7))
    data_train = spikeTrain[0:trainEndIdx]

    labels_train = dict()
    labels_train["isCorrect"] = spikeTrain_labels[0:trainEndIdx, 0]
    labels_train["isInSeq"] = spikeTrain_labels[0:trainEndIdx, 1]
    labels_train["trialPos"] = spikeTrain_labels[0:trainEndIdx, 2]
    labels_train["odor"] = spikeTrain_labels[0:trainEndIdx, 3]
    labels_train["trialID"] = trialID[0:trainEndIdx]

    data_test = spikeTrain[trainEndIdx:]

    labels_test = dict()
    labels_test["isCorrect"] = spikeTrain_labels[trainEndIdx:, 0]
    labels_test["isInSeq"] = spikeTrain_labels[trainEndIdx:, 1]
    labels_test["trialPos"] = spikeTrain_labels[trainEndIdx:, 2]
    labels_test["odor"] = spikeTrain_labels[trainEndIdx:, 3]
    labels_test["trialID"] = trialID[trainEndIdx:]

    dataset["inputTrain"] = data_train
    dataset["labelsTrain"] = labels_train

    dataset["inputTest"] = data_test
    dataset["labelsTest"] = labels_test

    dataset["times"] = times

    return dataset


# Get subset of data according to selected labels
def get_selected_data(labels_select, acc_measure, input_data, labels):
    idxs = np.ones(input_data.shape[0], dtype=bool)
    name = acc_measure

    odorNames = ["A", "B", "C", "D", "E"]
    trialNames = range(1, 6)

    for label in labels_select.keys():
        idxs_label = np.zeros(input_data.shape[0], dtype=bool)
        if label in ["odor", "trialPos"]:
            name = "%s_%s" % (name, label)

        for val in labels_select[label]:
            idxs_label = idxs_label | (labels[label] == val)

            if label == "odor":
                name = "%s%s" % (name, odorNames[val])
            elif label == "trialPos":
                name = "%s%s" % (name, trialNames[val])

        if label == "isInSeq":
            if labels_select[label][0] == 1:
                name = "%s_%s" % (name, "inSeq")
            else:
                name = "%s_%s" % (name, "outSeq")
        elif label == "isCorrect":
            if labels_select[label][0] == 1:
                name = "%s_%s" % (name, "correct")
            else:
                name = "%s_%s" % (name, "incorrect")

        idxs = idxs & idxs_label

    input_data_select = input_data[idxs]
    labels_select = labels[acc_measure][idxs]

    return idxs, input_data_select, labels_select, name
