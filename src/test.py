#!/usr/bin/env python3

"""
Created on April 12th 2018

@author: Chen
"""

import numpy as np
import argparse
from textwrap import dedent
from itertools import product
from keras.models import Sequential
from keras.layers import Dense, Embedding, Flatten
from keras import optimizers


# One hot encoding for bases
BASES = {'A': 0,  # 0001
         'C': 1,  # 0010
         'G': 2,  # 0100
         'T': 3,  # 1000
         'a': 0,
         'c': 1,
         'g': 2,
         't': 3}
VERSION = "v1.0.0"


def one_hot_encoding(seq, k_list, kmer_encoding):
    encoded_list = np.zeros(len(seq) * len(k_list) - sum(k_list) + len(k_list))
    seq = seq.upper()
    for k_i in range(len(k_list)):
        for i in range(len(seq) - k_list[k_i] + 1):
            kmer = seq[i: i + k_list[k_i]]
            encoded_list[i + len(seq) * k_i - sum(k_list[:k_i]) + k_i] = kmer_encoding[kmer]
    return encoded_list


def fasta_to_vectors(in_fasta, k_list):
    with open(in_fasta) as f:
        header_seq = f.readlines()

    seq = [header_seq[i * 2 + 1].strip() for i in range(int(len(header_seq)/2))]

    # Generate all unique kmers and their one hot encoding
    unique_kmers = sum(pow(4, k_list))
    all_kmers = []
    for k in k_list:
        all_kmers.extend(list(product(['A', 'T', 'C', 'G'], repeat=k)))
    all_kmers = [''.join(x) for x in all_kmers]
    kmer_encoding = dict(zip(all_kmers, range(unique_kmers)))

    seq_vector = [one_hot_encoding(x, k_list, kmer_encoding) for x in seq]
    return seq_vector


def create_model(l, weights=''):
    model = Sequential()

    nodes = 128
    k = np.array([4, 6, 8, 10])
    length = l * len(k) - sum(k) + len(k)
    model.add(Embedding(sum(pow(4, k)), nodes, input_length=length))
    model.add(Flatten())
    model.add(Dense(512, activation='relu'))
    model.add(Dense(64, activation='relu'))
    model.add(Dense(32, activation='relu'))
    model.add(Dense(3, activation='softmax'))

    adam = optimizers.Adam(lr=0.001)

    if weights:
        model.load_weights(weights)

    model.compile(loss='categorical_crossentropy',
                  optimizer=adam,
                  metrics=['accuracy'])
    return model


def main():
    parser = argparse.ArgumentParser(
        description=dedent('''
        Terminitor test
        -----------------------------------------------------------
        Given pre-trained models, predict whether a given sequence
        contains a poly(A) site at specific location or not
        '''),
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-v', '--version', action='version', version='Terminitor ' + VERSION)
    parser.add_argument('-t', '--test_file', help="Fasta file to be tested", required=True)
    parser.add_argument('-m', '--model', help="Pre-trained model file", required=True)
    parser.add_argument('-l', help="Length of input sequences", required=True, type=int)
    parser.add_argument('-o', '--output', help="Output probabilities", required=True)

    kmers = [4, 6, 8, 10]
    k = np.array(kmers)

    args = parser.parse_args()

    test_file = args.test_file
    model = args.model
    l = args.l
    out_f = args.output

    seq_vector = fasta_to_vectors(test_file, k)
    x_test = np.array(seq_vector)

    model = create_model(l, model)

    predict_prob = model.predict(x_test)
    out =  open(out_f, 'w')
    for i in range(len(x_test)):
        out.write(str(predict_prob[i][0]) + '\n')
    out.close()


if __name__ == "__main__":
    main()
