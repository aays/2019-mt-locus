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
        self.zstart1 = int(zstart1) # origin-zero
        self.end1 = int(end1)
        self.name2 = name2
        self.strand2 = strand2
        self.size2 = size2
        self.zstart2 = int(zstart2) # origin-zero, orientation dependent
        self.end2 = int(end2)
        self.identity = [int(num) for num in str(identity).split('/')]
        self.idPct = float(idPct[:len(idPct) - 1])
        self.coverage = [int(num) for num in str(coverage).split('/')]
        self.covPct = covPct[:len(covPct) - 1]

def parse_aln(filename):
    '''
    (str) -> list
    parses LASTZ alignment file.
    duplicates should have been removed using R script
    expects --format=general
    '''
    with open(filename) as f:
        aln_file = [aln(*line.split('\t')) 
                    for line in f.readlines() 
                    if not line.startswith('score')]
    return aln_file

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

def create_plus_dicts(filename, plus_length):
    '''
    (str, int) -> dict, dict

    input fasta should contain mt+ strains only

    creates two dictionaries:

    plus_seqs - sequences containing all Ns that
    will be written over with gametologs from
    the lastz alignment

    plus_refs - full mt+ sequences from fasta files, to
    be used to pull homologous regions from
    '''    
    plus_strains = [s.id for s in SeqIO.parse(filename, 'fasta')]
    plus_seqs = dict.fromkeys(plus_strains, '')
    plus_refs = dict.fromkeys(plus_strains, '')
    # N dictionary
    for strain in plus_seqs.keys():
        plus_seqs[strain] = ''.join(['N' for i in range(plus_length)])
    # ref dictionary
    for record in SeqIO.parse(filename, 'fasta'):
        plus_refs[record.id] = str(record.seq)
    return plus_seqs, plus_refs

def create_minus_dicts(filename, plus_length):
    '''
    (str, int) -> dict, dict, dict

    input fasta should contain mt- strains only

    creates three dictionaries:

    minus_seqs - 'minus' sequences containing all Ns -
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

def add_homologous_regions_plus(plus_seqs, aln_file, plus_refs, plus_length):
    '''
    (dict, aln, dict) -> dict
    takes in dictionaries from create_plus_dicts and
    iterates over alignment, pasting on homologous
    regions into N-only dictionary

    will also check that regions are not edited twice

    this function doesn't require a reverse complement
    dict since lastz always reports target coordinates
    in the 5' -> 3' orientation
    '''
    # make dummy dict to keep track of changes
    edit_check = dict.fromkeys(list(plus_seqs.keys()), '')
    for strain in edit_check.keys():
        edit_check[strain] = ''.join(['0' for i in range(plus_length)])
    added_count = 0

    # iterate through alignment
    for region in tqdm(aln_file):
        start, end = region.zstart1, region.end1
        for strain in plus_seqs.keys():
            left_chunk = plus_seqs[strain][0:start]
            added_chunk = plus_refs[strain][start:end]
            right_chunk = plus_seqs[strain][end:len(plus_seqs[strain])]
            plus_seqs[strain] = left_chunk + added_chunk + right_chunk
            added_count += len(added_chunk)

            # check for double edit
            try:
                current_region = [int(site) for site in list(edit_check[strain][start:end])]
                current_region = [site + 1 for site in current_region]
                assert 2 not in current_region # would indicate double edit
            except AssertionError:
                print('Error - one or more sites in the mt+ sequences')
                print('was edited twice. This indicates overlapping regions.')
                print('The offending region was', start, '-', end)
                sys.exit(1)
            else:
                left_edit_chunk = edit_check[strain][0:start]
                added_chunk = ''.join(['1' for i in range(end - start)])
                right_edit_chunk = edit_check[strain][end:len(edit_check[strain])]
                edit_check[strain] = left_edit_chunk + added_chunk + right_edit_chunk
    print(added_count / len(plus_seqs.keys()))
    return plus_seqs


def add_homologous_regions_minus(minus_seqs, aln_file, minus_refs, minus_rev_refs, plus_length):
    ''' 
    (dict, aln, dict, dict) -> dict
    takes in dictionaries from create_minus_dicts and
    iterates over alignment, pasting on homologous
    regions into N-only dictionary

    will also check that regions are not edited twice
    '''
    # make dummy dict to keep track of changes
    edit_check = dict.fromkeys(list(minus_seqs.keys()), '')
    for strain in edit_check.keys():
        edit_check[strain] = ''.join(['0' for i in range(plus_length)])
    added_count = 0

    # iterate through alignment
    for region in tqdm(aln_file):
        start, end = region.zstart1, region.end1
        start_minus, end_minus, orientation = region.zstart2, region.end2, region.strand2
        assert orientation in ['+', '-']
        for strain in minus_seqs.keys():
            if orientation == '+':
                left_chunk = minus_seqs[strain][0:start]
                added_chunk = minus_refs[strain][start_minus:end_minus] # homolog from minus ref
                right_chunk = minus_seqs[strain][end:len(minus_seqs[strain])]
                minus_seqs[strain] = left_chunk + added_chunk + right_chunk
                added_count += len(added_chunk)
            elif orientation == '-':
                left_chunk = minus_seqs[strain][0:start]
                added_chunk = minus_rev_refs[strain][start_minus:end_minus] # homolog from minus rev ref
                right_chunk = minus_seqs[strain][end:len(minus_seqs[strain])]
                minus_seqs[strain] = left_chunk + added_chunk + right_chunk
                added_count += len(added_chunk)
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
                left_edit_chunk = edit_check[strain][0:start]
                added_chunk = ''.join(['1' for i in range(end - start)])
                right_edit_chunk = edit_check[strain][end:len(edit_check[strain])]
                edit_check[strain] = left_edit_chunk + added_chunk + right_edit_chunk

    print(added_count / len(minus_seqs.keys()))
    return minus_seqs

def main():
    plus, minus, alignment, output = args()

    # parse aln file and get region length
    aligned_regions = parse_aln(alignment)
    mt_plus_length = get_mt_plus_length(plus)

    # create dicts
    plus_seqs, plus_refs = create_plus_dicts(plus, mt_plus_length)
    minus_seqs, minus_refs, minus_rev_refs = create_minus_dicts(minus, mt_plus_length)
    
    print('Extracting mt+ gametologs...')
    plus_seqs = add_homologous_regions_plus(plus_seqs,
        aligned_regions, plus_refs, mt_plus_length)
    print('Extracting mt- gametologs...')
    minus_seqs = add_homologous_regions_minus(minus_seqs, 
        aligned_regions, minus_refs, minus_rev_refs, mt_plus_length)

    print('Writing to file...')
    with open(output, 'w') as f:
        for strain in plus_seqs.keys():
            f.write('>' + strain + '\n')
            f.write(plus_seqs[strain] + '\n')
        for strain in minus_seqs.keys():
            f.write('>' + strain + '\n')
            f.write(minus_seqs[strain] + '\n')
    print('Done.')
    print('Hooray!')

if __name__ == '__main__':
    main()

        
