__author__ = 'Fule Liu'


import subprocess
import operator
from kmer import make_kmer_list


"""Prepare for PseKNC."""


class Protein():
    def __init__(self, seq_name, seq_desc, seq_content):
        self.seq_name = seq_name
        self.seq_desc = seq_desc
        self.seq_content = seq_content

    def get_seq_name(self):
        return self.seq_name

    def get_seq_desc(self):
        return self.seq_desc

    def get_seq_content(self):
        return self.seq_content

    def __str__(self):
        return "%s\t%s\n%s\n" % (self.seq_name, self.seq_desc, self.seq_content)


class RemoteProtein(Protein):
    def __init__(self, seq_name, seq_desc, seq_content, remote_acids, remote_profiles):
        super(RemoteProtein, self).__init__(seq_name, seq_desc, seq_content)
        self.remote_acids = remote_acids
        self.remote_profiles = remote_profiles

    def set_seq_content(self, content):
        self.seq_content = content

    def get_remote_acids(self):
        return self.remote_acids

    def get_remote_profiles(self):
        return self.remote_profiles

    def get_top_n_gram_acids(self, n):
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
        if n <= 0:
            raise ValueError("N is a integer larger than 0.")
        if n > len(self.remote_acids):
            raise ValueError("N is invalid, the n in top-n-gram cannot be larger than %d" % len(self.remote_acids))

        return ["".join(acids_n_gram[:n]) for acids_n_gram in zip(*[acids for acids in self.remote_acids])]

    def get_top_n_gram_profile1(self, n):
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
        if n <= 0:
            raise ValueError("N is a integer larger than 0.")
        if n > len(self.remote_profiles):
            raise ValueError("N is invalid, the n in top-n-gram cannot be larger than %d" % len(self.remote_profiles))

        return [sum(profiles_n_gram[:n]) for profiles_n_gram in zip(*[profiles for profiles in self.remote_profiles])]

    def __str__(self):
        print_remote_acids = "\n".join(self.remote_acids)
        print_remote_profiles = "\n".join([str(remote_profile) for remote_profile in self.remote_profiles])
        return "%s\t%s\n%s\n%s\n%s\n" % (self.seq_name, self.seq_desc, self.seq_content,
                                         print_remote_acids, print_remote_profiles)


def generate_profile(profile_jar_path, input_file, output_folder):
    """Use java -jar command to generate profile.

    Parameters
    ----------
    profile_jar_path: string
                      the Top-n-gram.jar path.
    input_file: string
                the fasta protein sequence file path.
    output_folder: string
                   the output profiles folder path.
    """
    cmd = "java -jar %s produce_frequency %s %s" % (profile_jar_path, input_file, output_folder)
    subprocess.Popen(cmd).wait()


def generate_top_n_gram_profile():
    pass


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
    for i in range(0, len_lines, 6):
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


def pseknc(remote_proteins, n, w, lamada, alphabet, theta_type=1):
    """This is a complete process in PseKNC."""
    kmers = [make_kmer_list(k, alphabet) for k in range(1, n+1)]
    vec = []

    for remote_protein in remote_proteins:
        # Get the normalized occurrence frequency in the protein sequence.
        f = []
        thetas = []
        for k in range(1, n+1):
            acid_fre = {}
            for kmer in kmers[k-1]:
                acid_fre[kmer] = 0

            remote_acids = remote_protein.get_top_n_gram_acids(k)
            for acid in remote_acids:
                acid_fre[acid] += 1
            fre_sum = float(sum(acid_fre.values()))
            sorted_acid_vals = sorted(acid_fre.items(), key=operator.itemgetter(0))
            # print(sorted_acid_vals)
            f.extend([e[1] / fre_sum for e in sorted_acid_vals])

            if k == n:
                thetas.extend(get_parallel_factor(k, lamada, remote_protein))

        # print(thetas)
        theta_sum = sum(thetas)
        denominator = n + w * theta_sum

        # temp_vec = [round(fre / denominator, 3) for fre in f]
        temp_vec = [fre / denominator for fre in f]
        for theta in thetas:
            # temp_vec.append(round(w * theta / denominator, 3))
            temp_vec.append(w * theta / denominator)

        vec.append(temp_vec)

    return vec


def get_parallel_factor(k, lamada, remote_protein):
    """Get the corresponding factor theta list."""
    thetas = []
    remote_profiles = remote_protein.get_top_n_gram_profile1(k)
    l = len(remote_profiles)

    for i in range(1, lamada + 1):
        temp_sum = 0.0
        for j in range(0, l - k - i + 1):
            temp_sum += (remote_profiles[j] - remote_profiles[j+i])**2
        thetas.append(temp_sum / (l - k - i + 1))

    return thetas


if __name__ == "__main__":
    remote_proteins = read_top_n_gram_file("Top-N2-gram2.txt", 2)

    for i, e in enumerate(remote_proteins):
        # print(e)
        # print(e.get_top_n_gram_acids(2))
        # print(e.get_top_n_gram_profile1(2))
        remote_proteins[i].set_seq_content("".join(e.get_top_n_gram_acids(2)))
        print(remote_proteins[i].get_seq_content())

    alphabet = "ACDEFGHIKLMNPQRSTVWY"

    w = 0.05
    lamada = 2
    res = pseknc(remote_proteins, 2, w, lamada, alphabet)
    for e in res:
        print(e)

    # with open('res.txt', 'w') as f:
    #     for e in res:
    #         print(e)
    #         print(sum(e))
    #         print(len(e))
    #         f.write(str(e))
    #         f.write('\n')

    # print(len(remote_proteins))

    # list_test = ["AAA", "BBB", "CCC"]
    # print("\n".join((list_test)))

    # list_test = [['123', '45'], ['67', '890']]
    # print("\n".join([str(e) for e in list_test]))

    # print(remote_proteins[0].get_seq_name())