#!/usr/bin/env python3

"""
Created on Mar 11th 2019

@author: Chen
"""

import pysam
import re
import pybedtools
import sys
from copy import deepcopy


def is_polya(seq, polya):
    if forward:
        if seq == 'A' * len(seq):
            return True
        elif seq.count('A') >= len(seq) / 2:
            return True
        else:
            return False
    else:
        if seq == 'T' * len(seq):
            return True
        elif seq.count('T') >= len(seq) / 2:
            return True
        else:
            return False


def rev_comp(seq):
    BASE = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    new_seq = seq[::-1]
    new_seq = ''.join([BASE[x] for x in new_seq])
    return new_seq


CHR = ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
       'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
       'chr20', 'chr21', 'chr22', 'chrX', 'chrY')


def main():
    trans = sys.argv[1]
    out_file = sys.argv[2]
    up_len = int(sys.argv[3])
    polya_len = int(sys.argv[4])

    out = open(out_file, 'w')

    with open(trans, 'r') as f:
        for line in f:
            if line[0] == '>':
                name = line.split()[0]
            else:
                transcript = line.strip()
                if transcript[0: polya_len] == 'T' * polya_len:
                    for i in range(len(transcript)):
                        if transcript[i] != 'T':
                            seq = rev_comp(transcript[i: i + up_len])
                            out.write(name + '\n' + seq + '\n')
                            break
                        else:
                            i += 1
                elif transcript[-polya_len:] == 'A' * polya_len:
                    for i in range(len(transcript) - 1, 0, -1):
                        if transcript[i] != 'A':
                            seq = transcript[i - up_len + 1: i + 1]
                            out.write(name + '\n' + seq + '\n')
                            break
                        else:
                            i += 1
                else:
                    continue
    out.close()

if __name__ == "__main__":
    main()

