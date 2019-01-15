'''
make_mt_only.py - provided fastas and lastz alignment, 
returns non-gametolog regions only from input fastas

takes in a fasta file (containing strains of one mating type)
as well as lastz output (--format=general) to return one fasta
file containing the same strains but with gametologous regions
masked with Ns

assumes that in lastz output, the mt+ was the target and 
the mt- the query
'''

import argparse
from tqdm import tqdm
from Bio import SeqIO
import re
import sys

def args():
    parser = argparse.ArgumentParser(description = 'Create fasta file with shared regions masked',
                                     usage = 'make_mt_only.py [options]')

    parser.add_argument('-f', '--fasta', required = True,
                        type = str, help = 'FASTA file containing strain-specific sequences.')
    parser.add_argument('-a', '--alignment', required = True,
                        type = str, help = 'LASTZ alignment output (--format=general).')
    parser.add_argument('-m', '--mt_allele', required = True,
                        type = str, help = 'Which mating type allele? [plus/minus]')
    parser.add_argument('-o', '--output', required = True,
                        type = str, help = 'File to write to.')

    args = parser.parse_args()

    return args.fasta, args.alignment, args.mt_allele, args.output

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
                    if not line.startswith('#score')]
    return aln_file

def get_sequence_length(filename):
    '''
    (str) -> int
    gets length of region from file.
    also checks that all are of the same length.

    equivalent to get_mt_plus_length from align_mt_fasta.py
    '''
    seq_length = 0
    counter = 0
    for record in SeqIO.parse(filename, 'fasta'):
        ind_length = len(str(record.seq))
        seq_length += len(str(record.seq))
        counter += 1
        assert seq_length > 0
    assert int(seq_length / counter) == int(ind_length) # make sure all seqs same length
    return ind_length

def get_non_shared_bases_plus(lastz_file, seq_length):
    '''
    (aln_file, str) -> list

    returns a list of positions that _do not_ have matches
    in the lastz alignment (ie are non-shared regions)

    specifically for plus allele (bc forward orientation only)
    '''
    intervals = [[region.zstart1, region.end1] 
                  for region in lastz_file]
    bases_covered = []

    for start, end in intervals:
        values = [num for num in range(start, end)]
        bases_covered.extend(values)

    all_bases = set([i for i in range(seq_length)]) # all possible values
    non_shared_bases = all_bases.difference(set(bases_covered))

    return sorted(list(set(non_shared_bases)))

def get_non_shared_bases_minus(lastz_file, seq_length):
    '''
    (aln_file, str) -> list

    returns a list of positions that _do not_ have matches
    in the lastz alignment (ie that are non-shared regions)

    all rev comp positions are converted to the positive orientation

    specifically for minus allele (bc both orientations need
    to be considered here)
    '''
    strands = ['+', '-']
    intervals = dict.fromkeys(strands, [])
    bases_covered = dict.fromkeys(strands, [])
    all_bases = dict.fromkeys(strands, set([i for i in range(seq_length)]))
    non_shared_bases = dict.fromkeys(strands)

    for key in intervals.keys():
        intervals[key] = [[region.zstart2, region.end2] 
                         for region in lastz_file if region.strand2 == key]
        for start, end in intervals[key]:
            values = [num for num in range(start, end)]
            bases_covered[key].extend(values)
        non_shared_bases[key] = sorted(list(
            all_bases[key].difference(set(bases_covered[key]))
            )
        )

    # convert rev orientation to fwd + combine
    non_shared_rev_fixed = set([seq_length - base - 1
                            for base in non_shared_bases['-']])
    non_shared_final = set(non_shared_bases['+']) \
                        .union(non_shared_rev_fixed)

    return non_shared_final # combined coordinates, all in fwd orientation

def create_mt_dicts(filename, seq_length):
    '''
    (str, int) -> dict, dict
    
    creates two dictionaries:

    strain_seqs - sequences containing all Ns that
    will be written over with gametologs from
    the lastz alignment

    strain_refs - full mt+ sequences from fasta files, to
    be used to pull homologous regions from
    '''    
    strainlist = [s.id for s in SeqIO.parse(filename, 'fasta')]
    strain_seqs = dict.fromkeys(strainlist, '')
    strain_refs = dict.fromkeys(strainlist, '')
    # N dictionary
    for strain in strain_seqs.keys():
        strain_seqs[strain] = ''.join(['N' for i in range(seq_length)])
    # ref dictionary
    for record in SeqIO.parse(filename, 'fasta'):
        strain_refs[record.id] = str(record.seq)
    return strain_seqs, strain_refs

def add_nonhomologous_regions(strain_seqs, non_shared_bases, strain_refs):
    '''
    (dict, list, dict) -> dict
    takes in dictionaries from create_mt_dicts and non-shared
    bases from strain-specific fxn -
    iterates over alignment, pasting on nonhomologous
    regions into N-only dictionary with strains as keys

    will also check that regions are not edited twice

    this function doesn't require a reverse complement
    dict since lastz always reports target (mt+) coordinates
    in the 5' -> 3' orientation and the mt- coordinates
    will have been converted to the positive orientation
    '''

    # make dummy dict to keep track of changes
    # works differently from the one in align_mt_fasta
    allele_length = len(strain_seqs[sorted(list(strain_seqs.keys()))[0]])
    edit_check = dict.fromkeys(list(strain_seqs.keys()), '')
    for strain in edit_check.keys(): # nested dict w/ positions as keys
        edit_check[strain] = dict.fromkeys([i for i in range(allele_length)], 0)

    # iterate through non-shared bases list
    for site in tqdm(non_shared_bases):
        for strain in strain_seqs.keys():
            left_chunk = strain_seqs[strain][0:site]
            added_site = strain_refs[strain][site]
            right_chunk = strain_seqs[strain][site + 1:len(strain_seqs[strain])]
            strain_seqs[strain] = left_chunk + added_site + right_chunk

        # check for double edit
            try:
                edit_check[strain][site] += 1
                assert edit_check[strain][site] < 2
            except AssertionError:
                print('Error - a site was edited twice')
                print('The offending site was', site)
                sys.exit(1)

    return strain_seqs

def main():
    fasta, alignment, mt_allele, output = args()

    # parse aln file and get region length
    aligned_regions = parse_aln(alignment)
    seq_length = get_sequence_length(fasta)

    if mt_allele == 'plus':
        print('mt+ allele selected.')
        print('Extracting non-shared mt+ sites...')
        non_shared_bases = get_non_shared_bases_plus(aligned_regions,
            seq_length)
    elif mt_allele == 'minus':
        print('mt- allele selected.')
        print('Extracting non-shared mt- sites')
        non_shared_bases = get_non_shared_bases_minus(aligned_regions,
            seq_length)

    strain_seqs, strain_refs = create_mt_dicts(fasta, seq_length)
    print('Creating gametolog-masked sequences...')
    strain_seqs = add_nonhomologous_regions(strain_seqs,
        non_shared_bases, strain_refs)

    print('Writing masked sequences to ', output, '.', sep = '')
    with open(output, 'w') as f:
        for strain in strain_seqs.keys():
            strain_name = re.search('CC[0-9]{4}', strain).group(0)
            f.write('>' + strain_name + '\n')
            f.write(strain_seqs[strain] + '\n')
    print('Done.')
    print('Hooray!')

if __name__ == '__main__':
    main()
        
