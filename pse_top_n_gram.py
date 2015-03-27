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
    def __init__(self, seq_name, seq_desc, seq_content, remote_acids, remote_profiles):
        super(RemoteProtein, self).__init__(seq_name, seq_desc, seq_content)
        self.remote_acids = remote_acids
        self.remote_profiles = remote_profiles

    def __str__(self):
        print_remote_acids = "\n".join(self.remote_acids)
        print_remote_profiles = "\n".join([str(remote_profile) for remote_profile in self.remote_profiles])
        return "%s\t%s\n%s\n%s\n%s\n" % (self.seq_name, self.seq_desc, self.seq_content,
                                         print_remote_acids, print_remote_profiles)

    def top_n_gram_acids(self, n):
        """Get top-n-gram acids.

        Parameter
        ---------
        n: int
           the n value in top-n-gram.

        Return
        ------
        top-n-gram acids: list
                          every elem represents the remote top-n-gram acid.
        """
        if n > len(self.remote_acids):
            print("N is invalid, the n in top-n-gram cannot be larger than %d" % len(self.remote_acids))
            return -1

        return ["".join(acids_n_gram[:n]) for acids_n_gram in zip(*[acids for acids in self.remote_acids])]

    def top_n_gram_profile1(self, n):
        """Get the sum of top-n-gram profile.

        Parameter
        ---------
        n: int
           the n value in top-n-gram.

        Return
        ------
        top-n-gram-profile: list
                            every elem represents the sum of the remote top-n-gram frequency value.
        """
        if n > len(self.remote_profiles):
            print("N is invalid, the n in top-n-gram cannot be larger than %d" % len(self.remote_profiles))
            return -1

        return [sum(profiles_n_gram) for profiles_n_gram in zip(*[profiles for profiles in self.remote_profiles])]


def read_top_n_gram_file(filename, n):
    """Get RemoteProtein from top-n-gram file.

    Parameter
    ---------
    filename: string
              the top-n-gram filename.
    n: int
       the n value in top-n-gram.

    Return
    ------
    remote_proteins: list
                     RemoteProtein objects.
    """
    with open(filename) as f_read:
        lines = f_read.readlines()

    remote_proteins = []
    len_lines = len(lines)
    step = 2 + 2*n
    for i in range(0, len_lines, step):
        seq_name = lines[i].strip()
        seq_desc = None
        seq_content = None
        remote_acids = []
        for j in range(i + 1, i + 2*n, 2):
            remote_acids.append(lines[j].strip())
        remote_profiles = []
        for j in range(i + 2, i + 2*n + 1, 2):
            remote_profiles.append([float(profile_val) for profile_val in lines[j].strip().split('\t')])

        remote_proteins.append(RemoteProtein(seq_name, seq_desc, seq_content, remote_acids, remote_profiles))

    return remote_proteins


def pseknc(input_data, k, w, lamada, phyche_list, alphabet, extra_index_file=None, all_prop=False, theta_type=1):
    """This is a complete process in PseKNC.

    :param k: int, the value of k-tuple.
    :param phyche_list: list, the input physicochemical properties list.
    :param extra_index_file: a file path includes the user-defined phyche_index.
    :param all_prop: bool, choose all physicochemical properties or not.
    """
    pass


if __name__ == "__main__":
    list_a = [[1, 2, 3, 4], [5, 6, 7, 8]]
    # print(len(list_a))
    # print(sum(list_a))
    #
    # protein = Protein(1, 2, 3)
    # print(protein)
    # protein2 = RemoteProtein(1, 2, 3, 4, 5)
    # print(protein2)
    #
    remote_proteins = read_top_n_gram_file("Top-N2-gram.txt", 2)
    for e in remote_proteins:
        print(e)
        print(e.top_n_gram_acids(1))
        print(e.top_n_gram_profile1(2))
    #
    # print(len(remote_proteins))

    # list_test = ["AAA", "BBB", "CCC"]
    # print("\n".join((list_test)))

    # list_test = [['123', '45'], ['67', '890']]
    # print("\n".join([str(e) for e in list_test]))


