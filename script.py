__author__ = 'Fule Liu'

import re


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


if __name__ == "__main__":
    delete_x("scop.1.53.1e-25.fasta", "scop.1.53.1e-25_delete_x.fasta")
    pass
