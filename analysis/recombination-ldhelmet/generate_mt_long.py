'''
generate_mt_long.py - generate mt locus in long format
'''

import argparse
from Bio import SeqIO
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(description = '',
                                     usage = 'script.py [options]')

    parser.add_argument('-m', '--mt_locus', required = True,
                        type = str, help = 'Aligned mt_locus LDhelmet')
    parser.add_argument('-p', '--plus', required = True,
                        type = str, help = 'Plus only LDhelmet')
    parser.add_argument('-a', '--alignment', required = True,
                        type = str, help = 'LASTZ bed file')
    parser.add_argument('-f', '--fasta', required = True,
                        type = str, help = 'Plus ref fasta')
    parser.add_argument('-o', '--outfile', required = True,
                        type = str, help = 'Outfile name')

    args = parser.parse_args()

    return args.mt_locus, args.plus, args.alignment, \
           args.fasta, args.outfile

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

def get_non_shared_bases(lastz_file, seq_length):
    '''
    (aln_file, str) -> list
    returns a list of positions that _do not_ have matches
    in the lastz alignment (ie are non-shared regions)
    specifically for plus allele (bc forward orientation only)
    '''
    intervals = [[region.zstart1 + 298299, region.end1 + 298299] 
                  for region in lastz_file]
    bases_covered = []

    for start, end in intervals:
        values = [num for num in range(start, end)]
        bases_covered.extend(values)

    all_bases = set([i + 298299 for i in range(seq_length)]) # all possible values
    non_shared_bases = all_bases.difference(set(bases_covered))

    return sorted(list(set(non_shared_bases)))

def ldhelmet_file(infile):
    '''
    returns ldhelmet file contents as dict
    '''
    out_dict = {}
    with open(infile, 'r') as f:
        for line in f:
            if line.startswith(('#', 'version', 'left')):
                continue
            else:
                line = line.split(' ')
                interval = [int(i) for i in line[0:2]]
                rho = float(line[2])
                for i in range(interval[0], interval[1]):
                    out_dict[i] = rho
    return out_dict
    

def write_long_mt(mt_locus, plus, intervals, outfile):
    mt_bases = range(298299, 943475) # ldhelmet is origin one
    mt_locus_dict = ldhelmet_file(mt_locus)
    plus_dict = ldhelmet_file(plus)
    
    with open(outfile, 'w') as f:
        f.write('position rho is_gametolog\n')
        for base in tqdm(mt_bases):
            line = str(base)
            if base in intervals: # ie non gametolog
                try:
                    line = line + ' ' + str(plus_dict[base]) + ' ' + '0'
                    f.write(line + '\n')
                    continue
                except KeyError:
                    line = line + ' NA 0'
                    f.write(line + '\n')
                    continue
            elif base not in intervals: # gametolog
                try:
                    line = line + ' ' + str(mt_locus_dict[base]) + ' ' + '1'
                    f.write(line + '\n')
                    continue
                except KeyError:
                    line = line + ' NA 1'
                    f.write(line + '\n')
                    continue

def main():
    mt_locus, plus, alignment, fasta, outfile = args()

    seq_length = get_sequence_length(fasta)
    aln_file = parse_aln(alignment)
    intervals = get_non_shared_bases(aln_file, seq_length)
    write_long_mt(mt_locus, plus, intervals, outfile)
    

if __name__ == '__main__':
    main()

        

