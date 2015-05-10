__author__ = 'Fule Liu'

import os
import subprocess

from pse_top_n_gram import read_top_n_gram_file, pseknc
import util

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"


def generate_vecs(top_n_gram_file, n, w, lamada):
    pos_train_remote_proteins = read_top_n_gram_file(top_n_gram_file, n)
    for i, e in enumerate(pos_train_remote_proteins):
        pos_train_remote_proteins[i].set_seq_content("".join(e.get_top_n_gram_acids(2)))
    res = pseknc(pos_train_remote_proteins, 2, w, lamada, ALPHABET)
    return res


def statistic(n, w, lamada, test_labels_file, predict_file, statistic_file):
    with open(test_labels_file) as f:
        lines = f.readlines()
    test_labels = {}
    for line in lines[1:]:
        line = line.rstrip().split()
        test_labels[line[0]] = line[1]

    # Get the predicted test labels.
    with open(predict_file) as f:
        lines = f.readlines()
    predict_score_label = {}
    for line in lines[12:]:
        # print(line)
        line = line.rstrip().split()
        predict_score_label[line[0]] = (line[2], test_labels[line[0]])

    # Write the statistic.
    roc_auc = util.roc_auc(predict_score_label)
    roc50_auc = util.roc50_auc(predict_score_label)
    tp, fp, tn, fn, tpr, fpr, precision = util.statistical(predict_score_label)

    with open(statistic_file, 'wb') as f:
        f.write(bytes("\t".join(["n:", str(n), '\n']), encoding="UTF-8"))
        f.write(bytes("\t".join(["w:", str(w), '\n']), encoding="UTF-8"))
        f.write(bytes("\t".join(["lamada:", str(lamada), '\n']), encoding="UTF-8"))
        f.write(bytes("\t".join(["ROC_AUC:", str(roc_auc), '\n']), encoding="UTF-8"))
        f.write(bytes("\t".join(["ROC50_AUC:", str(roc50_auc), '\n']), encoding="UTF-8"))
        f.write(bytes("\t".join(["TP:", str(tp), '\n']), encoding="UTF-8"))
        f.write(bytes("\t".join(["FP:", str(fp), '\n']), encoding="UTF-8"))
        f.write(bytes("\t".join(["TN:", str(tn), '\n']), encoding="UTF-8"))
        f.write(bytes("\t".join(["FN:", str(fn), '\n']), encoding="UTF-8"))
        f.write(bytes("\t".join(["TPR:", str(tpr), '\n']), encoding="UTF-8"))
        f.write(bytes("\t".join(["FPR:", str(fpr), '\n']), encoding="UTF-8"))
        f.write(bytes("\t".join(["Precision:", str(precision), '\n']), encoding="UTF-8"))


GIST_TRAIN_SVM = "gist/gist-train-svm.exe"
GIST_CLASSIFY = "gist/gist-classify.exe"

if __name__ == "__main__":
    listdir = os.listdir("data/dataset/")

    for cur_dir in listdir:
        tar_fold = "data/dataset/" + cur_dir + "/"
        pos_train = tar_fold + "pos_train"
        neg_train = tar_fold + "neg_train"
        pos_test = tar_fold + "pos_test"
        neg_test = tar_fold + "neg_test"
        tar_train = tar_fold + "train"
        tar_test = tar_fold + "test"
        tar_train_label = tar_fold + "train_label"
        tar_test_label = tar_fold + "test_label"

        # Generate pos_train vecs.
        n = 2
        w = 0.05
        lamada = 2
        pos_train_vecs = generate_vecs(pos_train, n, w, lamada)
        neg_train_vecs = generate_vecs(neg_train, n, w, lamada)
        pos_test_vecs = generate_vecs(pos_test, n, w, lamada)
        neg_test_vecs = generate_vecs(neg_test, n, w, lamada)
        train_vecs = pos_train_vecs + neg_train_vecs
        test_vecs = pos_test_vecs + neg_test_vecs
        train_labels = ['1'] * len(pos_train_vecs) + ['-1'] * len(neg_train_vecs)
        test_labels = ['1'] * len(pos_test_vecs) + ['-1'] * len(neg_test_vecs)

        util.write_gistsvm_vec(train_vecs, tar_train)
        util.write_gistsvm_vec(test_vecs, tar_test)
        util.write_gistsvm_class(train_labels, tar_train_label)
        util.write_gistsvm_class(test_labels, tar_test_label)

        # Train and predict.
        tar_weights = tar_fold + "weights"
        cmd = os.path.normcase("".join([GIST_TRAIN_SVM, " -train ", tar_train, " -class ", tar_train_label,
                                        " > ", tar_weights]))
        subprocess.Popen(cmd, shell=True).wait()
        tar_predict = tar_fold + "predict"
        cmd = os.path.normcase("".join([GIST_CLASSIFY, " -train ", tar_train, " -learned ", tar_weights,
                                        " -test ", tar_test, " > ", tar_predict]))
        subprocess.Popen(cmd, shell=True).wait()

        # Statistic.
        statistic_file = tar_fold + "statistic"
        statistic(n, w, lamada, tar_test_label, tar_predict, statistic_file)

        print("%s is ok.\n" % tar_fold)

    print("End!")

    # subprocess.Popen("gist/gist-train-svm.exe -train data/dataset/1.27.1.1/train -class data/dataset/1.27.1.1/train_label > data/dataset/1.27.1.1/weights").wait()
    # print("END.")