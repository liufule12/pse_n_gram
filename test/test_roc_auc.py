__author__ = 'Fule Liu'

import sys
sys.path.append("..")

import util


if __name__ == "__main__":
    # Get the test labels.
    with open("test.labels.txt") as f:
        lines = f.readlines()
    test_labels = {}
    for line in lines[1:]:
        line = line.rstrip().split()
        test_labels[line[0]] = line[1]

    # Get the predicted test labels.
    with open("test.predict") as f:
        lines = f.readlines()
    predict_score_label = {}
    for line in lines[12:]:
        # print(line)
        line = line.rstrip().split()
        predict_score_label[line[0]] = (line[2], test_labels[line[0]])

    print(predict_score_label)
    print(util.roc_auc(predict_score_label))
    print(util.statistical(predict_score_label))