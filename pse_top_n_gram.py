__author__ = 'Fule Liu'

import sys
import os
import pickle
from math import pow

import const
from util import frequency
from util import get_data
from util import check_args, read_k
from kmer import make_kmer_list


"""Prepare for PseKNC."""


class Protein():
    def __init__(self, seq_name, seq_desc, seq_content):
        self.seq_name = seq_name
        self.seq_desc = seq_desc
        self.seq_content = seq_content

    def __str__(self):
        return "%s\t%s\n%s\n" % (self.seq_name, self.seq_desc, self.seq_content)


class RemoteProtein(Protein):
    def __init__(self, seq_name, seq_desc, seq_content, remote_acids, remote_profile):
        super(RemoteProtein, self).__init__(seq_name, seq_desc, seq_content)
        self.remote_acids = remote_acids
        self.remote_profile = remote_profile

    def __str__(self):
        return "%s\t%s\n%s\n%s\n%s\n" % (self.seq_name, self.seq_desc, self.seq_content,
                                         self.remote_acids, self.remote_profile)

    def top_n_gram_acids(self, n):
        """Get top-n-gram acids."""
        if n > len(self.remote_acids):
            print("N is invalid, the n in top-n-gram cannot be larger than %d" % len(self.remote_acids))
            return -1

        return self.remote_acids[:n]

    def top_n_gram_profile1(self, n):
        """Get the sum of top-n-gram profile."""
        if n > len(self.remote_profile):
            print("N is invalid, the n in top-n-gram cannot be larger than %d" % len(self.remote_profile))
            return -1

        return sum(self.remote_profile)


def pseknc(input_data, k, w, lamada, phyche_list, alphabet, extra_index_file=None, all_prop=False, theta_type=1):
    """This is a complete process in PseKNC.

    :param k: int, the value of k-tuple.
    :param phyche_list: list, the input physicochemical properties list.
    :param extra_index_file: a file path includes the user-defined phyche_index.
    :param all_prop: bool, choose all physicochemical properties or not.
    """
    pass


if __name__ == "__main__":
    list_a = [1, 2, 3, 4]
    print(len(list_a))
    print(sum(list_a))

    protein = Protein(1, 2, 3)
    print(protein)
    protein2 = RemoteProtein(1, 2, 3, 4, 5)
    print(protein2)