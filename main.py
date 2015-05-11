__author__ = 'Fule Liu'

import os
import subprocess

from pse_top_n_gram import read_top_n_gram_file, pseknc
import util

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"


def generate_vecs(top_n_gram_file, n, w, lamada):
    """Generate pse_top_n_gram vecs."""
    pos_train_remote_proteins = read_top_n_gram_file(top_n_gram_file, n)
    for i, e in enumerate(pos_train_remote_proteins):
        pos_train_remote_proteins[i].set_seq_content("".join(e.get_top_n_gram_acids(n)))
    res = pseknc(pos_train_remote_proteins, n, w, lamada, ALPHABET)
    return res


def write_statistic(n, w, lamada, test_labels_file, predict_file, statistic_file):
    """Write the statistic result."""
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


def statistic_all(n, lamada, w, foldpath):
    """Statistic all statistic file result."""
    tp = 0
    fp = 0
    tn = 0
    fn = 0
    tpr = 0
    fpr = 0
    precision = 0
    roc_auc = 0
    roc50_auc = 0

    n_lamada_w = "_".join([str(n), str(lamada), str(w)])
    listdir = os.listdir(foldpath)
    for cur_dir in listdir:
        statistic_file = foldpath + cur_dir + "/statistic_" + n_lamada_w
        with open(statistic_file) as f:
            lines = f.readlines()
            roc_auc += float(lines[3].rstrip().split('\t')[1])
            roc50_auc += float(lines[4].rstrip().split('\t')[1])
            tp += float(lines[5].rstrip().split('\t')[1])
            fp += float(lines[6].rstrip().split('\t')[1])
            tn += float(lines[7].rstrip().split('\t')[1])
            fn += float(lines[8].rstrip().split('\t')[1])
            tpr += float(lines[9].rstrip().split('\t')[1])
            fpr += float(lines[10].rstrip().split('\t')[1])
            precision += float(lines[11].rstrip().split('\t')[1])


GIST_TRAIN_SVM = "gist/gist-train-svm.exe"
GIST_CLASSIFY = "gist/gist-classify.exe"

if __name__ == "__main__":
    listdir = os.listdir("data/dataset/")

    for cur_dir in listdir:
        n = 2
        for lamada in range(8, 10):
            w = 0.0
            while w < 1:
                w += 0.05

                # File path.
                n_lamada_w = "_".join([str(n), str(lamada), str(w)])
                tar_fold = "data/dataset/" + cur_dir + "/"
                pos_train = "".join([tar_fold, "pos_train"])
                neg_train = "".join([tar_fold, "neg_train"])
                pos_test = "".join([tar_fold, "pos_test"])
                neg_test = "".join([tar_fold, "neg_test"])
                tar_train = "".join([tar_fold, "train_", n_lamada_w])
                tar_test = "".join([tar_fold, "test_", n_lamada_w])
                tar_train_labels = "".join([tar_fold, "train_labels_", n_lamada_w])
                tar_test_labels = "".join([tar_fold, "test_labels_", n_lamada_w])
                tar_weights = "".join([tar_fold, "weights_", n_lamada_w])
                tar_predict = "".join([tar_fold, "predict", n_lamada_w])
                statistic_file = "".join([tar_fold, "statistic_", n_lamada_w])

                # Generate pos_train vecs.
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
                util.write_gistsvm_class(train_labels, tar_train_labels)
                util.write_gistsvm_class(test_labels, tar_test_labels)

                # Train and predict.
                cmd = os.path.normcase("".join([GIST_TRAIN_SVM, " -train ", tar_train, " -class ", tar_train_labels,
                                                " > ", tar_weights]))
                subprocess.Popen(cmd, shell=True).wait()
                cmd = os.path.normcase("".join([GIST_CLASSIFY, " -train ", tar_train, " -learned ", tar_weights,
                                                " -test ", tar_test, " > ", tar_predict]))
                subprocess.Popen(cmd, shell=True).wait()

                # Statistic.
                write_statistic(n, w, lamada, tar_test_labels, tar_predict, statistic_file)

                print("%s%s is ok.\n" % (tar_fold, n_lamada_w))

    print("End!")