__author__ = 'Fule Liu'

import re
import os


def delete_x(file_read, file_write):
    with open(file_read) as f_read:
        lines = f_read.readlines()
    with open(file_write, 'w') as f_write:
        for line in lines:
            if line[0] == '>':
                f_write.write(line)
                continue
            for e in line:
                if e != 'x' and e != 'X':
                    f_write.write(e)


def statistic_all(foldpath):
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

    listdir = os.listdir(foldpath)
    for cur_dir in listdir:
        if os.path.isfile(cur_dir):
            print(cur_dir)
            continue
        statistic_file = "/".join([foldpath, cur_dir, "statistic"])
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

    len_listdir = len(listdir)
    print(len_listdir)
    return tp/len_listdir, fp/len_listdir, tn/len_listdir, fn/len_listdir, tpr/len_listdir, fpr/len_listdir, \
           precision/len_listdir, roc_auc/len_listdir, roc50_auc/len_listdir

if __name__ == "__main__":
    # delete_x("scop.1.53.1e-25.fasta", "scop.1.53.1e-25_delete_x.fasta")
    print(statistic_all("data/dataset"))
    pass
