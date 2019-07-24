#!/usr/bin/env python3

"""
Created on Dec 11th, 2017

@author: Chen Yang

This script filters KLEAT results and compares with ground truth (simulated reads).
Generates metrics such as precision and sensitivity.

"""

import sys
import argparse


CHR = {"chr1": 1, "chr10": 2, "chr11": 3, "chr12": 4, "chr13": 5, "chr14": 6, "chr15": 7,
       "chr16": 8, "chr17": 9, "chr18": 10, "chr19": 11, "chr2": 12, "chr20": 13, "chr21": 14,
       "chr22": 15, "chr3": 16, "chr4": 17, "chr5": 18, "chr6": 19, "chr7": 20, "chr8": 21,
       "chr9": 22, "chrX": 23, "chrY": 24}


def comparison(infile, ref, tpm, fp_out, tp_out, fn_out):
    TP = 0
    TN = 0
    FP = 0
    FN = 0

    with open(infile, 'r') as f1, open(ref, 'r') as f2:
        kleat = f1.readline().split()
        simulated = f2.readline().split()
        while kleat and simulated:
            # print(TP, TN, FP, FN)
            # print(kleat, simulated)
            k_chr = kleat[0]
            s_chr = simulated[0]

            if k_chr not in CHR:
                kleat = f1.readline().split()
                continue
            elif s_chr not in CHR:
                simulated = f2.readline().split()
                continue

            if CHR[k_chr] > CHR[s_chr]:
                reads = float(simulated[4])
                if reads <= tpm:
                    TN += 1
                else:
                    FN += 1
                    fn_out.write('\t'.join(simulated) + '\n')
                simulated = f2.readline().split()
            elif CHR[k_chr] < CHR[s_chr]:
                FP += 1
                fp_out.write('\t'.join(kleat) + '\n')
                kleat = f1.readline().split()
            else:
                k_cs = int(kleat[1])
                s_cs = int(simulated[1])
                reads = float(simulated[4])

                if abs(k_cs - s_cs) <= 30 and kleat[5] == simulated[5]:
                    # print('\t'.join(kleat))
                    if reads > tpm:
                        TP += 1
                        tp_out.write('\t'.join(kleat) + '\t' + str(s_cs) + '\n')
                    kleat = f1.readline().split()
                    simulated = f2.readline().split()

                elif k_cs < s_cs:
                    FP += 1
                    fp_out.write('\t'.join(kleat) + '\n')
                    kleat = f1.readline().split()
                else:
                    if reads <= tpm:
                        TN += 1
                    else:
                        FN += 1
                        fn_out.write('\t'.join(simulated) + '\n')
                    simulated = f2.readline().split()
        while kleat:
            FP += 1
            fp_out.write('\t'.join(kleat) + '\n')
            kleat = f1.readline().split()
        while simulated:
            reads = float(simulated[4])
            if reads <= tpm:
                TN += 1
            else:
                FN += 1
                fn_out.write('\t'.join(simulated) + '\n')
            simulated = f2.readline().split()

    return TP, TN, FP, FN


def main(argv):
    # Parse input and output files
    parser = argparse.ArgumentParser()
    parser.add_argument("--i1", help="forward bed file")
    parser.add_argument("--i2", help="reverse bed file")
    parser.add_argument("--r1", help="forward reference bed file")
    parser.add_argument("--r2", help="reverse reference bed file")
    parser.add_argument('tpm', type=float, action="store")
    args = parser.parse_args(argv)

    infile1 = args.i1
    infile2 = args.i2
    ref1 = args.r1
    ref2 = args.r2
    tpm = args.tpm
    # print(tpm)

    if infile1 == '' or infile2 == '' or ref1 == '' or ref2 == '':
        print("Please specify the KLEAT bed file and reference bed file!")
        sys.exit(1)

    fp_out = open("FP", 'w')
    tp_out = open("TP", 'w')
    fn_out = open("FN", 'w')

    forward = comparison(infile1, ref1, tpm, fp_out, tp_out, fn_out)
    reverse = comparison(infile2, ref2, tpm, fp_out, tp_out, fn_out)

    TP = forward[0] + reverse[0]
    TN = forward[1] + reverse[1]
    FP = forward[2] + reverse[2]
    FN = forward[3] + reverse[3]

    sensitivity = TP / (TP + FN)
    #specificity = TN / (TN + FP)
    precision = TP / (TP + FP)
    #f1 = 2 * TP / (2 * TP + FP + FN)
    print("Sensitivity: " + str(sensitivity))
    # print("Specificity: " + str(specificity))
    print("Precision: " + str(precision))
    # print("F1: " + str(f1))
    print("TP: " + str(TP))
    print("TN: " + str(TN))
    print("FP: " + str(FP))
    print("FN: " + str(FN))

    fp_out.close()
    tp_out.close()
    fn_out.close()


if __name__ == "__main__":
    main(sys.argv[1:])
