'''
r2_calc.py - calculate r2 from long form sequences
generated by transpose_aligned_fasta.py
'''

import argparse
import csv
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(description = 'Calculate r2 stats',
                                     usage = 'python3.5 r2_calc.py [options]')

    parser.add_argument('-f', '--filename', required = True,
                        type = str, help = 'Input long-form file')
    parser.add_argument('-o', '--outfile', required = True,
                        type = str, help = 'File to write to')
    parser.add_argument('-w', '--windowsize', required = True,
                        type = str, help = 'Maximum dist b/w SNPs')

    args = parser.parse_args()

    return args.filename, args.outfile, int(args.windowsize)

    
def get_freqs(record1, record2):
    ''' (dict, dict) -> dict, dict
    returns two dicts:
    1. haplotype frequencies b/w two records (AB, ab, Ab, aB)
    2. allele frequencies (p1, p2, q1, q2)

    since there's no 'reference', the A and a alleles are
    arbitrarily assigned - this affects the sign of the D value,
    but D is squared when calculating r2 anyways

    if we have a non-gametolog SNP, ie
    A A T A N N N N (n = 8, first 4 is mt+, last 4 is mt-)

    and a SNP without variation in the mt+ only strains:
    C C C C G G C G

    then there are only two possible haps = AC and TC - in this
    case the fxn will return False

    this behaviour is checked for in usable_pair()
    '''
    record1_strains = [key for key in list(record1.keys()) 
                       if key != 'position' and record1[key] != 'N']
    record2_strains = [key for key in list(record1.keys()) 
                       if key != 'position' and record2[key] != 'N']
    usable_strains = set(record1_strains).intersection(set(record2_strains))

    try:
        A, a = list(set([record1[key] for key in usable_strains]))
        B, b = list(set([record2[key] for key in usable_strains]))
    except ValueError: # 'not enough values to unpack'
        return False

    haps_count = dict.fromkeys(['AB', 'ab', 'Ab', 'aB'], 0)
    freqs_count = dict.fromkeys(['p1', 'p2', 'q1', 'q2'], 0)
    total_calls = 0
    
    for strain in usable_strains:
        gt1 = record1[strain]
        gt2 = record2[strain]
        if gt1 == A and gt2 == B:
            haps_count['AB'] += 1
            freqs_count['p1'] += 1
            freqs_count['q1'] += 1
            total_calls += 1
        elif gt1 == a and gt2 == b:
            haps_count['ab'] += 1
            freqs_count['p2'] += 1
            freqs_count['q2'] += 1
            total_calls += 1
        elif gt1 == A and gt2 == b:
            haps_count['Ab'] += 1
            freqs_count['p1'] += 1
            freqs_count['q2'] += 1
            total_calls += 1
        elif gt1 == a and gt2 == B:
            haps_count['aB'] += 1
            freqs_count['p2'] += 1
            freqs_count['q1'] += 1
            total_calls += 1
    haps = dict.fromkeys(['AB', 'ab', 'Ab', 'aB'])
    freqs = dict.fromkeys(['p1', 'p2', 'q1', 'q2'])
    for key in haps:
        haps[key] = haps_count[key] / total_calls
    for key in freqs:
        freqs[key] = freqs_count[key] / total_calls

    return haps, freqs


def dcalc(record1, record2):
    ''' (dict, dict) -> float
    calculates Lewontin's D between pair of usable SNPs

    used as helper function for r2 calculation
    '''
    haps, freqs = get_freqs(record1, record2)
    try:
        lhs = haps['AB'] * haps['ab']
    except KeyError:
        lhs = 0
    try:
        rhs = haps['Ab'] * haps['aB']
    except KeyError:
        rhs = 0
    d = lhs - rhs
    return d

def r2calc(record1, record2):
    ''' (dict, dict) -> float
    calculates r2 between pair of usable SNPs
    '''
    haps, freqs = get_freqs(record1, record2)
    if 0 in freqs.values():
        r2 = 0
    else:
        d2 = dcalc(record1, record2) ** 2
        r2 = d2 / (freqs['p1'] * freqs['p2'] * \
                   freqs['q1'] * freqs['q2'])
    return round(r2, 5)

def usable_pair(record1, record2):
    '''
    checks whether the two records are usable
    for an LD calculation
    '''
    # usable SNP?
    if record2['is_snp'] == '0':
        return False

    # 'backwards' calculation/same SNP?
    record1['position'] = int(record1['position'])
    record2['position'] = int(record2['position'])
    if record1['position'] >= record2['position']:
        return False

    elif not get_freqs(record1, record2): # see get_freqs description
        return False
    
    else:
        return True

def ld_calc(infile, outfile, windowsize):
    '''(str, str, int) -> None
    calculates r2 between SNPs of specified windowsize
    and writes to an outfile
    '''
    with open(outfile, 'w') as f_out:
        f_out.write('snp1 snp2 r2\n')

        with open(infile, 'r') as f_in:
            ref_reader = csv.DictReader(f_in, delimiter = ' ')

            for record1 in tqdm(ref_reader):
                if record1['is_snp'] == '0':
                    continue
                with open(infile, 'r') as f_target:
                    target_reader = csv.DictReader(f_target, delimiter = ' ')

                    for record2 in target_reader:
                        usable = usable_pair(record1, record2)
                        distance = int(record1['position']) - int(record2['position'])
                        if distance > windowsize:
                            break
                        elif usable:
                            print(record1['position'], record2['position'])
                            ld_out = r2calc(record1, record2)
                            print(ld_out)
                            out_list = [record1['position'],
                                        record2['position'],
                                        ld_out]
                            out = ' '.join([str(i) for i
                                            in out_list])
                            f_out.write(out + '\n')
                        else:
                            continue

def main():
    infile, outfile, windowsize = args()
    ld_calc(infile, outfile, windowsize)

if __name__ == '__main__':
    main()

        

