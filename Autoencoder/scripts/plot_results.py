import numpy as np
import scipy.stats
import argparse

from utils import data_utils
import pylab
import os

import pickle


def rgb_to_hex(rgb_val):
    hexVal = '#%02x%02x%02x' % (rgb_val[0], rgb_val[1], rgb_val[2])
    return hexVal


def mean_sem_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, sem = np.mean(a), scipy.stats.sem(a)
    h = sem * scipy.stats.t.ppf((1 + confidence) / 2., n - 1)
    return m, m - sem, m + sem, m - h, m + h


def plot_encodings(data, encs, save_dir, dataset_name):
    # Define types of plots
    """
    labels_selects = [
                    dict(),
                    {"isInSeq":[1],"isCorrect":[1]},
                    {"odor":[0,1,2],"isInSeq":[1],"isCorrect":[1]},
                    {"isInSeq":[1],"isCorrect":[0]},
                    {"odor":[0,1,2],"isInSeq":[1],"isCorrect":[0]},

                    {"isCorrect":[1]},
                    {"isCorrect":[0]},

                    dict()
                   ]

    acc_measures = ["odor","odor","odor","odor","odor","isInSeq","isInSeq","isCorrect"]
    """
    labels_selects = [
        {"odor": [0, 1, 2, 3], "isInSeq": [1], "isCorrect": [1]},
        {"isCorrect": [1]},
        dict()
    ]

    acc_measures = ["odor", "isInSeq", "isCorrect"]

    # Define plot look
    # odorToColor = {0: 'black', 1: 'purple', 2:'blue', 3:'red', 4:'green'}
    odorToColor = {0: (44, 168, 224), 1: (154, 133, 122), 2: (9, 161, 74), 3: (128, 66, 151), 4: (241, 103, 36)}
    inSeqToColor = {0: (255, 0, 0), 1: (0, 0, 128)}
    correctToColor = {0: (255, 0, 0), 1: (0, 255, 0)}

    for key in odorToColor.keys():
        odorToColor[key] = rgb_to_hex(odorToColor[key])

    for key in inSeqToColor.keys():
        inSeqToColor[key] = rgb_to_hex(inSeqToColor[key])

    for key in correctToColor.keys():
        correctToColor[key] = rgb_to_hex(correctToColor[key])

    dotSize = 7
    alphaVal = 1.0

    # Iterate over each time window and plot encodings
    plotIdx = 1
    for enc_t, data_t in zip(encs, data):
        # Get labels
        enc = np.concatenate((enc_t["inputTrain"], enc_t["inputTest"]), axis=0)
        labels = dict()
        for labelName in data_t["labelsTrain"].keys():
            labels[labelName] = np.concatenate((data_t["labelsTrain"][labelName], data_t["labelsTest"][labelName]),
                                               axis=0)

        odorColors = np.array([odorToColor[val] for val in labels["odor"]])
        inSeqColors = np.array([inSeqToColor[val] for val in labels["isInSeq"]])
        correctColors = np.array([correctToColor[val] for val in labels["isCorrect"]])

        """
        inSeqIdxs = inSeqLabels == 1
        outSeqIdxs = inSeqLabels == 0
        odorABCInSeqIdxs = (inSeqLabels == 1) & (odorLabels <= 2)
        """

        startTime = min(data_t["times"])
        endTime = max(data_t["times"])
        # timeStr = "%.3fs to %.3fs" % (startTime,endTime)

        # Plot encoding for each selected subset and color coding
        for labels_select, acc_measure in zip(labels_selects, acc_measures):
            # Get selected subset
            idxs, encs_plot, labels_plot, name = data_utils.get_selected_data(labels_select, acc_measure, enc, labels)

            if acc_measure == "odor":
                plotColors = odorColors[idxs]
            elif acc_measure == "isInSeq":
                plotColors = inSeqColors[idxs]
            elif acc_measure == "isCorrect":
                plotColors = correctColors[idxs]
            else:
                raise ValueError("Unknown acc measure %s" % acc_measure)

            # Plot encodings
            fig, ax = pylab.subplots(1, 1)
            ax.scatter(encs_plot[:, 0], encs_plot[:, 1], color=plotColors, s=dotSize, alpha=alphaVal)

            # Plot centroids
            for lab in np.unique(plotColors):
                labIdxs = plotColors == lab
                plotEncs_mean = np.mean(encs_plot[labIdxs], axis=0)
                labColor = plotColors[labIdxs][0]

                ax.scatter(plotEncs_mean[0], plotEncs_mean[1], color=labColor, s=100 * dotSize, linewidths=5,
                           facecolors='white', alpha=0.5)
            # ax.set_title("%s %s" % (name,timeStr))

            # Save figure
            # saveDir_plt = "%s/%s" % (saveDir,name)
            # if not os.path.exists(saveDir_plt):
            #    os.mkdir(saveDir_plt)

            fig.savefig('%s/%s_%s_plot_%i_%ims_to_%ims.eps' % (save_dir, dataset_name, name, plotIdx, 1000 * startTime,
                                                               1000 * endTime))

            pylab.close(fig)

        plotIdx = plotIdx + 1


# Returns mat with dimensions numShuffleTrials x numTimepoints
def combine_shuff_accs(accs, key):
    acc = accs[0]
    shuff_acc_mat = np.stack([[accs[shuff_i][tp][key] for tp in range(len(acc))] for shuff_i in range(1, len(accs))],
                             axis=0)

    return shuff_acc_mat


def plot_acc(ax, fig, save_name, title, window_start_times, acc_vals, shuff_mat):
    # Plot accuracies
    color = 'k'
    ax.plot(window_start_times, acc_vals, color=color)

    # Plot randomly shuffled lable accuracies
    meanSEMCI = [None] * shuff_mat.shape[1]
    for tp in range(int(shuff_mat.shape[1])):
        meanSEMCI[tp] = mean_sem_confidence_interval(shuff_mat[:, tp], confidence=0.95)

    # Plot Mean
    color = 'w'
    shuff_mean = [x[0] for x in meanSEMCI]
    ax.plot(window_start_times, shuff_mean, color=color)

    # Plot confidence intervals
    color = 'm'
    CI_low = [x[3] for x in meanSEMCI]
    CI_high = [x[4] for x in meanSEMCI]
    ax.fill_between(window_start_times, CI_low, CI_high, color=color)

    # Plot standard error of the mean intervals
    color = 'k'
    SEMDiff_low = [x[1] for x in meanSEMCI]
    SEMDiff_high = [x[2] for x in meanSEMCI]
    ax.plot(window_start_times, SEMDiff_low, color=color, ls="--")
    ax.plot(window_start_times, SEMDiff_high, color=color, ls="--")

    # shuff_std = shuffMat.std(axis=0)
    # ax.fill_between(windowStartTimes,shuff_mean-shuff_std,shuff_mean+shuff_std,color=color)

    ax.legend()
    # ax.set_title(title)
    ax.set_xlabel('Time relative to nose-poke (s)')
    ax.set_ylabel('Accuracy')

    fig.savefig(save_name)


def plot_accs(accs, save_dir, dataset_name):
    plotNames = [x for x in accs[0][0].keys() if x != "startTime"]
    plotInfo = dict()
    for plotName in plotNames:
        plotInfo[plotName] = dict()
        plotInfo[plotName]["accVals"] = []

    windowStartTimes = []

    # Get times and accuracies
    for acc_t_idx, acc_t in enumerate(accs[0]):
        windowStartTimes.append(acc_t["startTime"])
        for plotName in plotNames:
            if not np.isnan(acc_t[plotName]):
                plotInfo[plotName]["accVals"].append(acc_t[plotName])
            else:
                plotInfo[plotName]["accVals"].append(0.0)

    # Plot accuracy
    for plotName in plotNames:
        saveName = '%s/%s_%s.eps' % (save_dir, dataset_name, plotName)
        fig, ax = pylab.subplots(1, 1)

        plot_acc(ax, fig, saveName, plotName, windowStartTimes, plotInfo[plotName]["accVals"],
                 combine_shuff_accs(accs, plotName))

        pylab.close(fig)

    print("")

    # Format plots
    # plotInfo["inSeqAcc"]["title"] = 'InSeq/OutSeq Accuracy'
    # plotInfo["odorAcc"]["title"] = 'Odor Accuracy'
    # plotInfo["odorBCDInSeqAcc"]["title"] = 'Odor (BCD,InSeq) Accuracy'
    # plotInfo["inSeqBCDAcc"]["title"] = 'InSeq/OutSeq (BCD) Accuracy'


def plot_results(save_dir_base, time_window, inp_len, num_shuffles, dataset_names):
    if not os.path.exists(save_dir_base):
        os.makedirs(save_dir_base)

    # Get data for training and testing
    agg_accs = None
    for dataset_name in dataset_names:
        print("Dataset: %s" % dataset_name)

        # Plot encodings
        data = data_utils.get_dataset(dataset_name, time_window, inp_len)
        encs = pickle.load(open("data/%s/encoding_%s_%s" % (dataset_name, time_window, inp_len), "rb"))

        plot_encodings(data, encs, save_dir_base, dataset_name)

        # Combine accuracies
        accs = pickle.load(open("data/%s/accs_%s_%s" % (dataset_name, time_window, inp_len), "rb"))
        if agg_accs is None:
            agg_accs = [None] * len(accs)

        for acc_idx, acc in enumerate(accs):
            if agg_accs[acc_idx] is None:
                agg_accs[acc_idx] = [None] * len(acc)
            for acc_t_idx, acc_t in enumerate(acc):
                if agg_accs[acc_idx][acc_t_idx] is None:
                    agg_accs[acc_idx][acc_t_idx] = dict()
                    for key in acc_t.keys():
                        agg_accs[acc_idx][acc_t_idx][key] = []
                for key in acc_t.keys():
                    if not np.isnan(acc_t[key]):
                        agg_accs[acc_idx][acc_t_idx][key].append(acc_t[key])

        # Plot accuracies
        plot_accs(accs, save_dir_base, dataset_name)

    # Plot aggregated accuracies
    for acc_idx, acc in enumerate(agg_accs):
        for acc_t_idx, acc_t in enumerate(acc):
            for key in acc_t.keys():
                agg_accs[acc_idx][acc_t_idx][key] = np.mean(agg_accs[acc_idx][acc_t_idx][key])

    plot_accs(agg_accs, save_dir_base, "all")


# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--results_dir', type=str, required=True, help="Location to save plots to")
parser.add_argument('--timewindow', type=float, required=True, help="Time window of plots")
parser.add_argument('--inpLen', type=int, required=True, help="Input length")
parser.add_argument('--num_shuffles', type=int, default=0, help="Number of random shuffles")
parser.add_argument('--datasets', type=str, required=True, help="Names of datasets separated by a comma")

args = parser.parse_args()

datasetNames = args.datasets.split(",")

plot_results(args.results_dir, args.timewindow, args.inpLen, args.num_shuffles, datasetNames)
