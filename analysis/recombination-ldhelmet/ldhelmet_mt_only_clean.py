'''
ldhelmet_mt_only_clean.py - 
remove gametologous region rhos from outfile

most functions from make_mt_only.py
'''

import argparse
from Bio import SeqIO
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(description = 'ldhelmet_mt_only_clean.py',
                                     usage = 'ldhelmet_mt_only_clean.py [options]')

    parser.add_argument('-i', '--filename', required = True,
                        type = str, help = 'LDhelmet infile')
    parser.add_argument('-b', '--bed', required = True,
                        type = str, help = 'Filtered LASTZ bed file')
    parser.add_argument('-a', '--allele', required = True,
                        type = str, help = 'mt allele [plus/minus]')
    parser.add_argument('-f', '--fasta', required = True,
                        type = str, help = 'FASTA file used as input to LDhelmet')
    parser.add_argument('-o', '--outfile', required = True,
                        type = str, help = 'File to write to')

    args = parser.parse_args()

    return args.filename, args.bed, args.allele, args.fasta, args.outfile

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

    return sorted(list(non_shared_final)) # combined coordinates, all in fwd orientation

def in_interval(line, intervals):
    '''(str, list) -> bool
    checks interval against a list of _non_ gametolog positions
    returns True if the interval is a subset of the larger position list
    '''
    line = line.split(' ')
    start, end = [int(i) for i in line[0:2]]
    current_interval = set(list(range(start - 1, end))) # LDhelmet is 1-based
    if not isinstance(intervals, set):
        intervals = set(intervals)
    return current_interval.issubset(intervals) # returns True if subset

def correct_line(line, offset):
    ''' (list) -> str
    returns modified line as a str
    '''
    line = line.split(' ')
    line[0] = int(line[0]) + offset
    line[1] = int(line[1]) + offset
    line = [str(i) for i in line]

    return ' '.join(line)

def filter_ldhelmet(infile, intervals, allele, outfile):
    with open(outfile, 'w') as out:
        if allele == 'plus':
            offset = 298298
            plus = True
        elif allele == 'minus':
            offset = 0
            plus = False

        # write header
        colnames = ['left_snp', 'right_snp', 'mean',
                    'p0.025', 'p0.500', 'p0.975\n']
        out.write(' '.join(colnames))

        with open(infile, 'r') as f:
            for line in tqdm(f):
                if line.startswith(('#', 'version')):
                    continue
                else:
                    if in_interval(line, intervals) and plus:
                        line = correct_line(line, offset)
                        out.write(line)
                    elif in_interval(line, intervals) and not plus:
                        out.write(line)
                    elif not in_interval(line, intervals):
                        continue

def main():
    filename, bed, allele, fasta, outfile = args()
    seq_length = get_sequence_length(fasta)

    aln_file = parse_aln(bed)

    if allele == 'plus':
        intervals = get_non_shared_bases_plus(aln_file, seq_length)
    elif allele == 'minus':
        intervals = get_non_shared_bases_minus(aln_file, seq_length)
    filter_ldhelmet(filename, intervals, allele, outfile)

if __name__ == '__main__':
    main()

