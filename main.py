__author__ = 'Fule Liu'


import kmer
from pse_top_n_gram import read_top_n_gram_file, pseknc
from util import write_gistsvm_vec, write_gistsvm_class

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"


def generate_vecs(top_n_gram_file, n, w, lamada):
    pos_train_remote_proteins = read_top_n_gram_file(top_n_gram_file, n)
    for i, e in enumerate(pos_train_remote_proteins):
        pos_train_remote_proteins[i].set_seq_content("".join(e.get_top_n_gram_acids(2)))
    res = pseknc(pos_train_remote_proteins, 2, w, lamada, ALPHABET)
    return res


if __name__ == "__main__":
    tar_fold = "data/dataset/7.3.5.2/"
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
    train_labels = ['1'] * len(pos_train_vecs) + ['-1'] * len(neg_train_vecs)
    test_vecs = pos_test_vecs + neg_test_vecs
    test_labels = ['1'] * len(pos_test_vecs) + ['-1'] * len(neg_test_vecs)

    write_gistsvm_vec(train_vecs, tar_train)
    write_gistsvm_vec(test_vecs, tar_test)
    write_gistsvm_class(train_labels, tar_train_label)
    write_gistsvm_class(test_labels, tar_test_label)

    print("End!")
    pass