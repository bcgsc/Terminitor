#!/usr/bin/env python3

"""
Created on Tue Aug 7th 2018

@author: Chen
"""

import sys
import numpy as np
import random
import argparse
from textwrap import dedent
from time import strftime
from itertools import product
from keras.models import Sequential
from keras.layers import Dense, Embedding, Flatten
from keras import optimizers
from sklearn.model_selection import KFold
from keras.callbacks import EarlyStopping, ModelCheckpoint


# One hot encoding for bases
BASES = {'A': 0,  # 0001
         'C': 1,  # 0010
         'G': 2,  # 0100
         'T': 3,  # 1000
         'a': 0,
         'c': 1,
         'g': 2,
         't': 3}
VERSION = 'v1.0.0'


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


def create_model(l, k, weights=''):
    model = Sequential()
    length = l * len(k) - sum(k) + len(k)
    model.add(Embedding(sum(pow(4, k)), 128, input_length=length))
    model.add(Flatten())

    model.add(Dense(512, activation='relu'))
    model.add(Dense(64, activation='relu'))
    model.add(Dense(32, activation='relu'))
    model.add(Dense(3, activation='softmax'))

    adam = optimizers.Adam(lr=0.001)

    if weights:
        model.load_weights(weights)
        print("Created model and loaded weights from file: ", weights)
    else:
        print(model.summary())
        print("adam: 0.001", )
    model.compile(loss='categorical_crossentropy',
                  optimizer=adam,
                  metrics=['accuracy'])

    return model


def main():
    parser = argparse.ArgumentParser(
        description=dedent('''
        Terminitor training
        -----------------------------------------------------------
        Given three fasta files with different labels: poly(A) CS, 
        non-poly(A) CS, and non-CS, train the model
        '''),
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-v', '--version', action='version', version='Terminitor ' + VERSION)
    parser.add_argument('-polya', help="Poly(A) CS, fasta file", required=True)
    parser.add_argument('-cs', help="Non-poly(A) CS, fasta file", required=True)
    parser.add_argument('-non', help="Non-CS, fasta file", required=True)
    parser.add_argument('-model', help="File name of trained model", required=True)
    parser.add_argument('-l', help="Length of input sequences", required=True, type=int)

    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Start reading fasta\n")
    sys.stdout.flush()

    args = parser.parse_args()

    lab1_fasta = args.polya
    lab2_fasta = args.cs
    lab3_fasta = args.non
    weight_file = args.model
    l = args.l

    k = np.array([4, 6, 8, 10])

    seq_vector_lab1 = np.array(fasta_to_vectors(lab1_fasta, k))
    seq_vector_lab2 = np.array(fasta_to_vectors(lab2_fasta, k))
    seq_vector_lab3 = np.array(fasta_to_vectors(lab3_fasta, k))
    lab_vector_lab1 = np.tile([1, 0, 0], (len(seq_vector_lab1), 1))
    lab_vector_lab2 = np.tile([0, 1, 0], (len(seq_vector_lab2), 1))
    lab_vector_lab3 = np.tile([0, 0, 1], (len(seq_vector_lab3), 1))

    # Build the whole model
    random.seed(123)

    size1 = len(seq_vector_lab1)
    size2 = len(seq_vector_lab2)
    size3 = len(seq_vector_lab3)

    train_i1 = random.sample(range(size1), int(0.8 * size1))
    train_i2 = random.sample(range(size2), int(0.8 * size2))
    train_i3 = random.sample(range(size3), int(0.8 * size3))
       
    val_i1 = [i for i in range(size1) if i not in train_i1]
    val_i2 = [i for i in range(size2) if i not in train_i2]
    val_i3 = [i for i in range(size3) if i not in train_i3]

    x_train = np.concatenate((seq_vector_lab1[train_i1], seq_vector_lab2[train_i2], seq_vector_lab3[train_i3]))
    y_train = np.concatenate((lab_vector_lab1[train_i1], lab_vector_lab2[train_i2], lab_vector_lab3[train_i3]))
    x_val = np.concatenate((seq_vector_lab1[val_i1], seq_vector_lab2[val_i2], seq_vector_lab3[val_i3]))
    y_val = np.concatenate((lab_vector_lab1[val_i1], lab_vector_lab2[val_i2], lab_vector_lab3[val_i3]))

    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Start training\n")
    sys.stdout.flush()

    model = create_model(l, k)
    early_stop = EarlyStopping(monitor='val_loss', patience=10)
    filepath = weight_file + "_best_weights.hdf5"
    checkpoint = ModelCheckpoint(filepath, monitor='val_loss', verbose=0, save_best_only=True)

    model.fit(x_train, y_train, batch_size=64, verbose=2, epochs=100,
              validation_data=(x_val, y_val),
              callbacks=[checkpoint, early_stop])
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Finished!\n")


if __name__ == "__main__":
    main()
