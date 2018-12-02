'''
align_mt_fasta_maf.py - create aligned fasta from maf file
'''

from tqdm import tqdm
import argparse
from Bio.Seq import Seq
from Bio import SeqIO

def args():
    parser = argparse.ArgumentParser(description = 'Create fasta file containing mt+ and mt- sequences',
                                     usage = 'align_mt_fasta_maf.py [options]')

    parser.add_argument('-p', '--plus', required = True,
                        type = str, help = 'FASTA file containing mt+ sequences.')
    parser.add_argument('-m', '--minus', required = True,
                        type = str, help = 'FASTA file containing mt- sequences.')
    parser.add_argument('-a', '--alignment', required = True,
                        type = str, help = 'LASTZ alignment output (--format=maf).')
    parser.add_argument('-o', '--output', required = True,
                        type = str, help = 'File to write to.')

    args = parser.parse_args()

    return [args.plus, args.minus, args.alignment, args.output]

class aln(object):
    '''
    Quick class to parse maf lastz alignment output.
    Expects --format=maf formatted output from command-line lastz
    '''
    def __init__(self, score, name1, zstart1, end1, strand1,
            size1, seq1, name2, zstart2, end2, strand2,
            size2, seq2):
        self.score = score
        self.name1 = name1
        self.zstart1 = int(zstart1) # origin-zero
        self.end1 = int(end1)
        self.strand1 = strand1
        self.size1 = size1
        self.seq1 = SeqIO.SeqRecord(Seq(seq1))
        self.name2 = name2
        self.zstart2 = int(zstart2) # origin-zero, orientation dependent
        self.end2 = int(end2)
        self.strand2 = strand2
        self.size2 = size2
        self.seq2 = SeqIO.SeqRecord(Seq(seq2))

class aln_file(object):
    def __init__(self, filename):
        with open(filename, 'r') as f:
            d = [line for line in f.readlines()
                 if not line.startswith('#')]

        alignments = []
        chunks = [d[i:i+4]
                  for i in range(0, len(d), 4)]

        for curr_aln in chunks:
            score = curr_aln[0].split('=')[1]
            for i in [1,2]: # sequences
                curr_aln[i] = [element for element 
                               in curr_aln[i].split(' ')
                               if element] # not empty

            name1, zstart1, end1, strand1, \
                size1, seq1 = curr_aln[1].split(' ')[1:]
            name2, zstart2, end2, strand2, \
                size2, seq2 = curr_aln[2].split(' ')[1:]

            aln_obj = aln(score, name1, zstart1, end1,
                strand1, size1, seq1, name2, zstart2,
                end2, strand2, size2, seq2)

            alignments.append(aln_obj)

        return alignments

# use biopython's mutable seqs to work with strain specific seqs!
# from Bio.Seq import MutableSeq
# myseq.tomutable()


