import argparse
import pickle
import csv
import numpy as np

import os

from utils import data_utils


def convert_results(time_window, inp_len, out_dir, dataset_names):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    for datasetName in dataset_names:
        print("Dataset: %s" % datasetName)

        outDir_ds = "%s/%s" % (out_dir, datasetName)
        if not os.path.exists(outDir_ds):
            os.mkdir(outDir_ds)

        #  Write encodings csv
        encs = pickle.load(open("data/%s/encoding_%s_%s" % (datasetName, time_window, inp_len), "rb"))
        data = data_utils.get_dataset(datasetName, time_window, inp_len)
        with open("%s/encoding_%s_%s.csv" % (outDir_ds, time_window, inp_len), mode='w') as csv_file:
            csv_writer = csv.writer(csv_file, delimiter=',')
            csv_writer.writerow(
                ['StartTime', 'CenterTime', 'TrialID', 'Trial', 'Odor', 'IsInSeq', 'IsCorrect', 'Dimension1',
                 'Dimension2'])

            for enc_t, data_t in zip(encs, data):
                startTime = min(enc_t["times"])
                startTime_str = "%.3f" % startTime
                centerTime_str = "%.3f" % (startTime + time_window / 2.0)
                enc = np.concatenate((enc_t["inputTrain"], enc_t["inputTest"]), axis=0)

                labels = dict()
                for labelName in data_t["labelsTrain"].keys():
                    labels[labelName] = np.concatenate(
                        (data_t["labelsTrain"][labelName], data_t["labelsTest"][labelName]), axis=0)

                for pointIdx, point in enumerate(enc):
                    csv_writer.writerow(
                        [startTime_str, centerTime_str, labels['trialID'][pointIdx], labels['trialPos'][pointIdx],
                         labels['odor'][pointIdx], labels['isInSeq'][pointIdx], labels['isCorrect'][pointIdx], point[0],
                         point[1]])

        # Write accuracies csv
        accs = pickle.load(open("data/%s/accs_%s_%s" % (datasetName, time_window, inp_len), "rb"))

        with open("%s/accs_%s_%s.csv" % (outDir_ds, time_window, inp_len), mode='w') as csv_file:
            csv_writer = csv.writer(csv_file, delimiter=',')
            csv_writer.writerow(['StartTime', 'odor', 'isInSeq', 'isCorrect', 'shuffleIdx'])

            for acc_idx, acc in enumerate(accs):
                for acc_t in acc:
                    csv_writer.writerow(["%.3f" % (acc_t['startTime']), acc_t['odor_inSeq_odorBCD_correct'],
                                         acc_t['isInSeq_odorBCD_correct'], acc_t['isCorrect'], acc_idx])


def main():
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--timewindow', type=float, required=True, help="Time window of plots")
    parser.add_argument('--inpLen', type=int, required=True, help="Input length")
    parser.add_argument('--datasets', type=str, required=True, help="Names of datasets separated by a comma")
    parser.add_argument('--out_dir', type=str, required=True, help="Output directory")

    args = parser.parse_args()

    datasetNames = args.datasets.split(",")

    convert_results(args.timewindow, args.inpLen, args.out_dir, datasetNames)


if __name__ == "__main__":
    main()
