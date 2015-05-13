__author__ = 'Fule Liu'

import os


def data_split(top_n_gram_file, membership_file, out_fold):
    with open(top_n_gram_file) as f:
        lines_top_n_gram = f.readlines()
    with open(membership_file) as f:
        lines_membership = f.readlines()

    # Create family dataset fold.
    families = lines_membership[0].rstrip().split()[1:]
    for family in families:
        fold_name = "/".join([out_fold, family])
        if not os.path.exists(fold_name):
            os.mkdir(fold_name)
        else:
            # print("Exist %s" % fold_name)
            pass

    # Construct dataset.
    for i, line in enumerate(lines_membership[1:]):
        line = line.rstrip().split()
        seq_name = line[0]
        for index, val in enumerate(line[1:]):
            print(seq_name)

            # According to all_seq_fam_membersship.txt construct dataset.
            if val == "1":
                filename = "/".join([out_fold, families[index], "pos_train"])
            elif val == "2":
                filename = "/".join([out_fold, families[index], "neg_train"])
            elif val == "3":
                filename = "/".join([out_fold, families[index], "pos_test"])
            elif val == "4":
                filename = "/".join([out_fold, families[index], "neg_test"])
            elif val == "0":
                continue
            else:
                print("Error val: %s" % val)
                print(seq_name)
                continue
            with open(filename, 'ab') as f:
                start_index = (i - 1) * 6
                end_index = i * 6
                for l in lines_top_n_gram[start_index: end_index]:
                    f.write(bytes(l, 'UTF-8'))

    print("End.")


if __name__ == "__main__":
    top_n_gram_file = "Top-N2-gram.txt"
    membership_file = "all_seq_fam_membership.txt"
    out_fold = "dataset"
    data_split(top_n_gram_file, membership_file, out_fold)