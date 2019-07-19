#!/usr/bin/env python3

"""
Created on April 12th 2018

@author: Chen
"""

import sys
import numpy as np
from time import strftime
from itertools import product
from keras.models import Sequential
from keras.layers import Dense, Embedding, Flatten, Dropout
from keras import optimizers
from keras.wrappers.scikit_learn import KerasClassifier
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
'''
BASES = {'A': [0, 0, 0, 1],  # 0001
         'C': [0, 0, 1, 0],  # 0010
         'G': [0, 1, 0, 0],  # 0100
         'T': [1, 0, 0, 0],  # 1000
         'a': [0, 0, 0, 1],
         'c': [0, 0, 1, 0],
         'g': [0, 1, 0, 0],
         't': [1, 0, 0, 0]}
'''


'''
def one_hot_encoding(seq, k=4):
    # Generate all unique kmers and their one hot encoding
    unique_kmers = pow(4, k)
    all_kmers = list(product(['A', 'T', 'C', 'G'], repeat=k))
    all_kmers = [''.join(x) for x in all_kmers]
    kmer_encoding = dict(zip(all_kmers, range(unique_kmers)))

    encoded_list = np.zeros(len(seq) - k + 1)
    seq = seq.upper()
    for i in range(len(seq) - k + 1):
        kmer = seq[i: i+k]
        encoded_list[i] = kmer_encoding[kmer]
    return encoded_list


def fasta_to_vectors(in_fasta, k=4):
    with open(in_fasta) as f:
        header_seq = f.readlines()

    seq = [header_seq[i * 2 + 1].strip() for i in range(int(len(header_seq)/2))]
    seq_vector = [one_hot_encoding(x, k) for x in seq]
    return seq_vector
'''


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

    '''
    model.add(Embedding(pow(4, 4), 48, input_length=197))
    # model.add(Dropout(0.2))
    # if rnn == "LSTM":
        # model.add(LSTM(nodes, input_shape=(200, 4)))
    # elif rnn == "GRU":
        # model.add(GRU(nodes, input_shape=(200, 4), activation='relu'))
    # model.add(Flatten())
    model.add(Dense(32, activation='relu'))
    # model.add(Dropout(0.5))
    model.add(Dense(32, activation='relu'))
    # model.add(Dropout(0.5))
    model.add(Dense(3, activation='softmax'))

    sgd = optimizers.SGD(lr=0.01, decay=0, momentum=0.9, nesterov=True)
    adagrad = optimizers.Adagrad(lr=0.01)
    rmsprop = optimizers.rmsprop(lr=0.001)
    adam = optimizers.Adam(lr=0.001)

    if weights:
        model.load_weights(weights)
        # print("Created model and loaded weights from file: ", weights)

    model.compile(loss='binary_crossentropy',
                  optimizer=sgd,
                  metrics=['accuracy'])
    # print(model.summary())
    '''

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
    test_fasta = sys.argv[1]
    filepath = sys.argv[2]
    l = int(sys.argv[3])

    kmers = [4, 6, 8, 10]
    k = np.array(kmers)
    seq_vector = fasta_to_vectors(test_fasta, k)
    x_test = np.array(seq_vector)

    model = create_model(l, filepath)

    predict_prob = model.predict(x_test)
    for i in range(len(x_test)):
        print(predict_prob[i][0])


if __name__ == "__main__":
    main()
