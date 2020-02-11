#!/usr/bin/env python

"""
Created on Apr 12th, 2018

@author: Chen Yang

This script combines features in bed file within certain distance.

"""

import sys
import getopt
import numpy as np
import operator


def cluster_helper(cs, d):
    coord_list = []

    if len(cs) == 0:
        return coord_list
    else:
        if len(set([x[0] for x in cs.values()])) == 1:
            coord = int(np.median(list(cs.keys())))
        else:
            coord = max(cs.items(), key=operator.itemgetter(1))[0]
        left_cs = {k: v for (k, v) in cs.items() if k < coord - d}
        right_cs = {k: v for (k, v) in cs.items() if k > coord + d}
        coord_list = cluster_helper(left_cs, d) + [coord] + cluster_helper(right_cs, d)
    return coord_list


def clustering(cs, d):
    '''
    :param cs: a dictionary, k = cleavage site, v = (prob, transcript)
    :return clusters: a dictionary with all representative cleavage site(s)
    '''

    coord_list = cluster_helper(cs, d)
    clusters = {}
    for coord in coord_list:
        if coord not in cs:
            prob = cs[list(cs.keys())[0]][0]
            trans = cs[list(cs.keys())[0]][1]
        else:
            prob = cs[coord][0]
            trans = cs[coord][1]
        clusters[coord] = [prob, trans]

    return clusters


def main():
    d = 0
    out_prefix = "out"

    # Parse options and parameters
    if len(sys.argv) < 2:
        print("Please specify a bed file!")
        sys.exit(1)
    else:
        try:
            opts, args = getopt.getopt(sys.argv[1:], "d:o:")
        except getopt.GetoptError:
            sys.exit(1)
        for opt, arg in opts:
            if opt == "-d":
                d = int(arg)
            elif opt == "-o":
                out_prefix = arg
            else:
                sys.exit(1)

    bed_file = sys.argv[-1]
    
    out = open(out_prefix + ".cluster.bed", 'w')
    # Merge features within distance d
    with open(bed_file, 'r') as f:
        last_line = f.readline()
        info = last_line.strip().split()
        last_chr = info[0]
        last_strand = info[5]
        cs = {int(info[1]): (float(info[4]), info[3])}
 
        for line in f:
            info = line.strip().split()
            if info[0] != last_chr or int(info[1]) >= max(cs.keys()) + d or info[5] != last_strand:
                clusters = clustering(cs, d)
                for coord, v in clusters.items():
                    prob, trans = v[0], v[1]
                    out_line = last_chr + '\t' + str(coord) + '\t' + str(coord + 1) + \
                               '\t' + trans + '\t' + str(prob) + '\t' + last_strand + '\n'
                    out.write(out_line)
                last_chr = info[0]
                last_strand = info[5]
                cs = {int(info[1]): (float(info[4]), info[3])}
            else:
                if int(info[1]) not in cs:
                    cs[int(info[1])] = (float(info[4]), info[3])
                elif cs[int(info[1])][0] < float(info[4]):
                    cs[int(info[1])] = (float(info[4]), info[3])
                
        if len(set([x[0] for x in cs.values()])) == 1:
            coord = int(np.median(list(cs.keys())))
        else:
            coord = max(cs.items(), key=operator.itemgetter(1))[0]
        if coord not in cs:
            prob = cs[list(cs.keys())[0]][0]
            trans = cs[list(cs.keys())[0]][1]
        else:
            prob = cs[coord][0]
            trans = cs[coord][1]

        out_line = last_chr + '\t' + str(coord) + '\t' + str(coord + 1) + \
                   '\t' + trans + '\t' + str(prob) + '\t' + last_strand + '\n'
        out.write(out_line)

    out.close()            


if __name__ == "__main__":
    main()
