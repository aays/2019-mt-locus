'''
align_mt_fasta.py - aligner that outputs a fasta file

takes in two fasta files (containing mt+ strains and mt- 
strains respectively) as well as lastz output (--format==general) 
to return one fasta file containing all mt+ records and mt-
regions homologous to the mt+
'''

import argparse
from tqdm import tqdm
from Bio import SeqIO
import sys

def args():
    parser = argparse.ArgumentParser(description = 'Create fasta file containing mt+ and mt- sequences',
                                     usage = 'align_mt_fasta.py [options]')

    parser.add_argument('-p', '--plus', required = True,
                        type = str, help = 'FASTA file containing mt+ sequences.')
    parser.add_argument('-m', '--minus', required = True,
                        type = str, help = 'FASTA file containing mt- sequences.')
    parser.add_argument('-a', '--alignment', required = True,
                        type = str, help = 'LASTZ alignment output (--format=general).')
    parser.add_argument('-o', '--output', required = True,
                        type = str, help = 'File to write to.')

    args = parser.parse_args()

    return [args.plus, args.minus, args.alignment, args.output]

class aln(object):
    '''
    Quick class to parse lastz alignment output.
    Expects --format=general formatted output from command-line lastz
    '''
    def __init__(self, score, name1, strand1, size1, zstart1, 
        end1, name2, strand2, size2, zstart2, end2, identity, 
        idPct, coverage, covPct):
        self.score = score
        self.name1 = name1
        self.strand1 = strand1
        self.size1 = size1
        self.zstart1 = zstart1 # origin-zero
        self.end1 = end1
        self.name2 = name2
        self.strand2 = strand2
        self.size2 = size2
        self.zstart2 = zstart2 # origin-zero, orientation dependent
        self.end2 = end2
        self.identity = [int(num) for num in str(identity).split('/')]
        self.idPct = float(idPct[:len(idPct) - 1])
        self.coverage = [int(num) for num in str(coverage).split('/')]
        self.covPct = float(covPct[:len(covPct) - 1])

def parse_aln(filename):
    '''
    (str) -> list
    parses LASTZ alignment file.
    duplicates should have been removed using R script
    expects --format=general
    '''
    with open(filename) as f:
        aln_file = [aln(line) for line in f.readlines() if not line.startswith('#')]
    return aln_file

def parse_mt_plus(filename, output):
    '''
    (str, str) -> None
    Writes mt plus sequences to final output file
    '''
    with open(output, 'w') as f:
        seqs = SeqIO.parse(filename, 'fasta')
        SeqIO.write(seqs, f, 'fasta') # write records to file
    return None

def get_mt_plus_length(filename):
    '''
    (str) -> int
    gets length of mt plus region from file.
    also checks that all are of the same length.
    '''
    plus_length = 0
    counter = 0
    for record in SeqIO.parse(filename, 'fasta'):
        ind_length = len(str(record.seq))
        plus_length += len(str(record.seq))
        counter += 1
        assert plus_length > 0
    assert int(plus_length / counter) == int(ind_length) # make sure all seqs same length
    return ind_length

def create_minus_dicts(filename, plus_length):
    '''
    (str, int) -> dict, dict, dict
    creates three dictionaries:

    minus_seqs - 'minus' sequences containing all Ns
    these sequences will be 'written over' with mt+
    homologous regions from the lastz alignment

    minus_refs - full mt- sequences from fasta files,
    to be used to pull homologous regions from 

    minus_rev_refs - same as minus_refs, but reversed
    '''
    minus_strains = [s.id for s in SeqIO.parse(filename, 'fasta')]
    minus_seqs = dict.fromkeys(minus_strains, '')
    minus_refs = dict.fromkeys(minus_strains, '')
    minus_rev_refs = dict.fromkeys(minus_strains, '')
    # N dictionary
    for strain in minus_seqs.keys():
        minus_seqs[strain] = ''.join(['N' for i in range(plus_length)])
    # ref dictionaries
    for record in SeqIO.parse(filename, 'fasta'):
        minus_refs[record.id] = str(record.seq)
        minus_rev_refs[record.id] = str(record.reverse_complement().seq)
    return minus_seqs, minus_refs, minus_rev_refs


def add_homologous_regions(minus_seqs, aln_file, minus_ref_dict,
    minus_ref_rev_dict):
    ''' (dict, aln, dict, dict) -> dict
    takes in dictionaries from create_minus_dicts and
    iterates over alignment, pasting on homologous
    regions into N-only dictionary

    will also check that regions are not edited twice
    '''
    # make dummy dict to keep track of changes
    plus_length = len(minus_seqs[sorted(list(minus_seqs.keys()))[0]])
    edit_check = dict.fromkeys(list(minus_seqs.keys()), '')
    for strain in edit_check.keys():
        edit_check[strain] = ''.join(['0' for i in range(plus_length)])

    # iterate through alignment
    for region in tqdm(aln_file):
        start, end, orientation = region.zstart2, region.end2, region.strand2
        assert orientation in ['+', '-']
        for strain in minus_dict.keys():
            if orientation == '+':
                minus_seqs[strain][start:end] = minus_refs[strain][start:end]
            elif orientation == '-':
                minus_seqs[strain][start:end] = minus_rev_refs[strain][start:end]
            # check for double edit
            try:
                current_region = [int(site) for site in list(edit_check[strain][start:end])]
                current_region = [site + 1 for site in current_region]
                assert 2 not in current_region # would indicate double edit
            except AssertionError:
                print('Error - one or more sites in the mt- sequences')
                print('was edited twice. This indicates overlapping regions.')
                print('The offending region was', start, '-', end, 'w/ orientation', orientation)
                sys.exit(1)
            else:
                edit_check[strain][start:end] = ''.join(['1' for i in range(end - start)])

    return minus_seqs

def main(plus, minus, alignment, output):
    aligned_regions = parse_aln(alignment)
    mt_plus_length = get_mt_plus_length(plus)
    minus_seqs, minus_refs, minus_rev_refs = create_minus_dicts(minus, mt_plus_length)

    # write initial output file to disk
    parse_mt_plus(plus, output) 

    minus_seqs = add_homologous_regions(minus_seqs, 
        aligned_regions, minus_refs, minus_rev_refs)

    with open(output, 'a') as f:
        for strain in minus_seqs.keys():
            f.write('>' + strain + '\n')
            f.write(minus_seqs[strain] + '\n')

if __name__ == '__main__':
    arguments = args()
    main(*arguments)

        
