#!/usr/bin/env python3

'''
Created on Sept 11th 2018
Major change on Oct 17th, adopting Design 9

@author: Chen
'''

import pysam
import pybedtools
import argparse
from HTSeq import GFF_Reader, parse_GFF_attribute_string
from textwrap import dedent
from copy import deepcopy


VERSION = 'v1.0.0'


def is_polya(seq, forward):
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


CHR = ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
       'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
       'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY')


def main():
    parser = argparse.ArgumentParser(
        description=dedent('''
        Terminitor pipeline extract candidate sequence
        -----------------------------------------------------------
        This script is the last step of Terminitor RNA-seq pipeline.
        It extracts candidate sequence from alignment bam file for
        testing by the pre-trained neural network model.
        '''),
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-v', '--version', action='version', version='Terminitor ' + VERSION)
    parser.add_argument('-t', '--annot_trans', help='Transcript annotation file, GTF format. This file contains '
                                                    'features whose type is transcript, can be downloaded from the '
                                                    'auxiliary files or created by users', required=True)
    parser.add_argument('-a', '--annot_all', help='Annotation file, GTF format. For example annotations from Ensembl '
                                                  'ftp site', required=True)
    parser.add_argument('-m', '--aln', help='The alignment file from assembled transcript contigs to reference genome in '
                                            'bam format.', required=True)
    parser.add_argument('-g', '--genome', help='Reference genome assembly in Fasta format. Can be downloaded from '
                                               'Ensembl ftp site', required=True)
    parser.add_argument('-o', help='Output file, fasta format containing candidate sequences to be tested',
                        required=True)
    parser.add_argument('-u', '--up_len', help='Upstream sequence length', type=int, default=100)
    parser.add_argument('-d', '--down_len', help='Downstream sequence length', type=int, default=100)

    args = parser.parse_args()

    annot_trans = args.annot_trans
    annot_all = args.annot_all
    aln = args.aln
    assem = args.genome
    out_file = args.o
    up_len = args.up_len
    down_len = args.down_len

    '''
    # Read in genome
    chrome_len = {}

    with open(assem, 'r') as f:
        for line in f:
            if line[0] == '>':
                info = line.split()
                name = info[0][1:]
                length = int(info[2].split(':')[4])
                chrome_len[name] = length
    '''

    # GTF file is 1-based inclusive, bed file is 0-based half open
    trans_coord = {}
    gtf_features = GFF_Reader(annot_all)
    for feature in gtf_features:
        # Sometimes, 'chr' is not in the chromesome name, add it back
        chrom = feature.iv.chrom if 'chr' in feature.iv.chrom else  'chr' + feature.iv.chrom
        if chrom not in CHR:
            continue
        if chrom not in trans_coord:
            trans_coord[chrom] = {}
        if feature.type == 'transcript':
            trans = feature.attr['transcript_id']
            trans_coord[chrom][trans] = {'exon': [], 'CDS': [], 'strand': feature.iv.strand}
        elif feature.type in ['exon', 'CDS'] :
            trans = feature.attr['transcript_id']
            coord = [feature.iv.start, feature.iv.end]
            trans_coord[chrom][trans][feature.type].append(coord)

    for chrom, transs in trans_coord.items():
        for trans, info in transs.items():
            info['exon'].sort()
            info['CDS'].sort()

            strand = info['strand']
            info['utr3'] = []
            info['utr5'] = []
            info['introns'] = []

            # Introns
            start = info['exon'][0][1]
            exon_p = 1
            while exon_p < len(info['exon']):
                info['introns'].append([start, info['exon'][exon_p][0]])
                start = info['exon'][exon_p][1]
                exon_p += 1

            # UTRs
            if info['CDS']:
                exon_p = 0
                CDS_p = 0
                while exon_p < len(info['exon']) and CDS_p < len(info['CDS']):
                    if info['exon'][exon_p][1] <= info['CDS'][CDS_p][0]:
                        if strand == '+':
                            info['utr5'].append(list(info['exon'][exon_p]))
                        else:
                            info['utr3'].append(list(info['exon'][exon_p]))
                        exon_p += 1
                        continue

                    if info['exon'][exon_p][0] < info['CDS'][CDS_p][0]:
                        if strand == '+':
                            info['utr5'].append([info['exon'][exon_p][0], info['CDS'][CDS_p][0]])
                        else:
                            info['utr3'].append([info['exon'][exon_p][0], info['CDS'][CDS_p][0]])

                    if info['CDS'][CDS_p][1] < info['exon'][exon_p][1]:
                        if strand == '+':
                            info['utr3'].append([info['CDS'][CDS_p][1], info['exon'][exon_p][1]])
                        else:
                            info['utr5'].append([info['CDS'][CDS_p][1], info['exon'][exon_p][1]])
                    CDS_p += 1
                    exon_p += 1
                while exon_p < len(info['exon']):
                    if strand == '+':
                        info['utr3'].append(list(info['exon'][exon_p]))
                    else:
                        info['utr5'].append(list(info['exon'][exon_p]))
                    exon_p += 1
                if info['utr3']:
                    if strand == '+' and info['utr3'][0][1] - info['utr3'][0][0] <= 3:
                        del info['utr3'][0]
                    elif strand == '-' and info['utr3'][-1][1] - info['utr3'][-1][0] <= 3:
                        del info['utr3'][-1]
                    else:
                        if strand == '+':
                            info['utr3'][0][0] += 3
                        else:
                            info['utr3'][-1][1] -= 3
            else:
                info['CDS'] = deepcopy(info['exon'])

    # Process alignment
    ensembl = pybedtools.BedTool(annot_trans)
    alignment = pybedtools.BedTool(aln)
    intersect = alignment.intersect(ensembl, bed=True, wo=True, split=True)

    seq_dict = {}

    for info in intersect:
        feature = parse_GFF_attribute_string(info[20][:-1] + '\n')
        trans = feature['transcript_id']
        chrom = info[0] if 'chr' in info[0] else 'chr' + info[0]
        if chrom not in CHR:
            continue
        
        strand = info[5]
        if strand == '+':
            name = info[3] + '_' + info[0] + '_' + str(int(info[2]) - 1) + '_' + 'F'
            dis2annot = abs(int(info[2]) - trans_coord[chrom][trans]['exon'][-1][1])
            if trans_coord[chrom][trans]['utr3']:
                utr3 = True if int(info[2]) - trans_coord[chrom][trans]['utr3'][0][0] > 0 else False
            else:
                utr3 = False
        else:
            name = info[3] + '_' + info[0] + '_' + info[1] + '_' + 'R'
            dis2annot = abs(int(info[1]) - trans_coord[chrom][trans]['exon'][0][0])
            if trans_coord[chrom][trans]['utr3']:
                utr3 = True if trans_coord[chrom][trans]['utr3'][-1][1] - int(info[1]) > 0 else False
            else:
                utr3 = False
        if name not in seq_dict:
            seq_dict[name] = {'trans': trans, 'dis2annot': dis2annot, 'utr3': utr3}
        else:
            if seq_dict[name]['utr3'] and utr3 and seq_dict[name]['dis2annot'] > dis2annot:
                seq_dict[name] = {'trans': trans, 'dis2annot': dis2annot, 'utr3': utr3}
            elif seq_dict[name]['utr3'] and (not utr3):
                continue
            elif (not seq_dict[name]['utr3']) and utr3:
                seq_dict[name] = {'trans': trans, 'dis2annot': dis2annot, 'utr3': utr3}
            elif (not seq_dict[name]['utr3']) and (not utr3) and seq_dict[name]['dis2annot'] > dis2annot:
                seq_dict[name] = {'trans': trans, 'dis2annot': dis2annot, 'utr3': utr3}

    samfile = pysam.AlignmentFile(aln, 'rb')
    for read in samfile.fetch(until_eof=True):
        cigar_string = read.cigartuples
        qseq = read.query_sequence
        if read.flag != 0 and read.flag != 16:
            continue
        
        chrom = read.reference_name if 'chr' in read.reference_name else 'chr' + read.reference_name
        if chrom not in CHR:
            continue

        if read.flag == 0:
            strand = '+'
            cs = read.reference_end - 1
            name = read.query_name + '_' + chrom + '_' + str(cs) + '_' + 'F'
        elif read.flag == 16:
            strand = '-'
            cs = read.reference_start
            name = read.query_name + '_' + chrom + '_' + str(cs) + '_' + 'R'

        if strand == '+':
            if cigar_string[-1][0] == 4 or cigar_string[-1][0] == 5:
                clipped = cigar_string[-1][1]
                clipped_seq = qseq[-clipped:]
                if not is_polya(clipped_seq, True):
                    continue
                u_seq = qseq[-clipped - up_len: -clipped]
            elif cigar_string[-1][0] == 0 and name in seq_dict and seq_dict[name]['utr3']:
                u_seq = qseq[-up_len:]
            else:
                continue
        
        else:
            if cigar_string[0][0] == 4 or cigar_string[0][0] == 5:
                clipped = cigar_string[0][1]
                clipped_seq = qseq[:clipped]
                if not is_polya(clipped_seq, False):
                    continue
                u_seq = qseq[clipped: clipped + up_len]
            elif cigar_string[0][0] == 0 and name in seq_dict and seq_dict[name]['utr3']:
                u_seq = qseq[:up_len]
            else:
                continue

        if name not in seq_dict:
            seq_dict[name] = {}
        seq_dict[name]['U'] = u_seq
        
    bed_string = []
    for k, v in seq_dict.items():
        info = k.split('_')
        if 'U' not in v:
            continue
        if info[-1] == 'F':
            bed_string.append(info[1] + '\t' + str(int(info[2]) + 1) + '\t' + str(int(info[2]) + down_len + 1) +
                              '\t' + k + '\t-\t+')
        else:
            bed_string.append(info[1] + '\t' + str(int(info[2]) - down_len) + '\t' + str(int(info[2])) + '\t' + k +
                              '\t-\t-')

    down_bed = pybedtools.BedTool('\n'.join(bed_string), from_string=True)

    down_seqs = down_bed.sequence(fi=assem, name=True, split=True)

    out = open(out_file, 'w')

    with open(down_seqs.seqfn) as f:
        for line in f:
            name = line[1:-1].split(':')[0]
            seq = next(f).strip()
            seq_dict[name]['D'] = seq
    nonredundant = {}    
    for name, v in seq_dict.items():
        direction = name.split('_')[-1]
        if 'U' not in v or len(v['U']) < up_len:
            continue
        if direction == 'F':
            seq = v['U'] + v['D']
        elif direction == 'R':
            seq = v['D'] + v['U']
            seq = rev_comp(seq)
        if 'N' in seq:
            continue
        if seq not in nonredundant:
            out.write('>' + name + '\n' + seq + '\n')
            nonredundant[seq] = 0

    out.close()


if __name__ == '__main__':
    main()

