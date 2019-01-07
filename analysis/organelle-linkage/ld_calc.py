'''
ld_calc.py - calculate LD stats between two regions in a single VCF
'''

import argparse
import vcf
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(description = 'ld_calc.py',
                                     usage = 'python3.5 ld_calc.py [options]')

    parser.add_argument('-v', '--vcf_file', required = True,
                        type = str, help = 'VCF file to calculate LD from.')
    parser.add_argument('-r', '--regions', required = True,
                        type = str, nargs = '+', help = 'Two inter-chr regions to calculate LD between.')
    parser.add_argument('-o', '--outfile', required = True,
                        type = str, help = 'File to write to.')

    args = parser.parse_args()

    return args.vcf_file, args.regions, args.outfile


def snppuller(vcf_file, chrom = None, start = None, end = None):
    vcfin = vcf.Reader(filename = vcf_file, compressed = True)

    def is_diallelic(record):
        counts = list(set([len(record.REF), len(record.ALT), len(record.ALT[0])]))
        if len(counts) == 1 and counts[0] == 1:
            return True
        else:
            return False

    def is_nonsingleton(record):
        if isinstance(record.INFO['AN'], list):
            count = record.INFO['AN'][0] - record.INFO['AC'][0]
        else:
            count = record.INFO['AN'] - record.INFO['AC'][0]
        if count != 1 and record.INFO['AC'][0] != 1:
            return True
        else:
            return False

    def is_variant(record):
        if record.INFO['AC'][0] != 0 and record.INFO['AF'][0] != 1.0:
            return True
        else:
            return False

    def usable_record(record):
        if is_diallelic(record) and is_nonsingleton(record) and is_variant(record):
            return True
        else:
            return False

    args = [item for item in [chrom, start, end] if item]
    for record in vcfin.fetch(*args):
        if usable_record(record):
            yield record
        else:
            continue

def straingetter(record1, record2, GQ_threshold = 30, DP_threshold = 20):
    rec1set = set([s.sample for s in record1.samples if s['GT'] != '.'])
    rec2set = set([s.sample for s in record2.samples if s['GT'] != '.'])
    strainlist = list(rec1set.intersection(rec2set))
    # filter for GQ
    garbage = []
    for strain in strainlist:
        gt1 = record1.genotype(strain)
        gt2 = record2.genotype(strain)
        if gt1['GQ'] < GQ_threshold or gt2['GQ'] < GQ_threshold:
            garbage.append(strain)
        elif gt1['DP'] < DP_threshold or gt2['DP'] < DP_threshold:
           garbage.append(strain)
    if garbage:
        for trash_strain in garbage:
            strainlist.remove(trash_strain)
    return strainlist

def get_freqs(record1, record2):
    strainlist = straingetter(record1, record2)
    haplist = []
    p_count, q_count, total_calls = 0, 0, 0
    haps = dict.fromkeys(['AB', 'ab', 'Ab', 'aB'], 0)
    for strain in strainlist:
        gt1 = record1.genotype(strain)['GT']
        gt2 = record2.genotype(strain)['GT']
        if '.' in [gt1, gt2]:
            continue
        elif gt1 == '0' and gt2 == '0': # AB
            haps['AB'] += 1
            p_count += 1
            q_count += 1
            total_calls += 1
        elif gt1 == '1' and gt2 == '1': # ab
            haps['ab'] += 1
            total_calls += 1
        elif gt1 == '1' and gt2 == '0': # aB
            haps['aB'] += 1
            q_count += 1
            total_calls += 1
        elif gt1 == '0' and gt2 == '1': # Ab
            haps['Ab'] += 1
            p_count += 1
            total_calls += 1
    values = dict.fromkeys(['p1', 'p2', 'q1', 'q2'], 0.0)
    if not total_calls:
        values['p1'], values['q1'], values['p2'], values['q2'] = 0, 0, 0, 0
        return values, haps 
    else:
        values['p1'] = p_count / total_calls
        values['q1'] = q_count / total_calls
        values['p2'] = 1 - values['p1']
        values['q2'] = 1 - values['q1']
        for key in haps:
            haps[key] = haps[key] / total_calls

        haplist = [list(key) for key in haps if haps[key] != 0]
        first = list(set([hap[0] for hap in haplist]))
        second = list(set([hap[1] for hap in haplist]))
        if len(first) == 1 or len(second) == 1:
            return False # one SNP is invariant
        else:
            return values, haps


def ld_calc(values, haps):
    '''(dict, dict) -> float, float, float
    calculates D, D', and r2 before returning all as a tuple
    '''
    # D
    try:
        LHS = haps['AB'] * haps['ab']
    except KeyError: # either hap missing
        LHS = 0
    try:
        RHS = haps['Ab'] * haps['aB']
    except KeyError:
        RHS = 0
    d = LHS - RHS

    # Lewontin's D'
    if d >= 0:
        dmax = min(values['p1'] * values['q2'], values['p2'] * values['q1'])
        if dmax == 0:
            dprime = 0
        else:
            dprime = d / dmax
    elif d < 0:
        dmin = max(-1 * values['p1'] * values['q1'], -1 * values['p2'] * values['q2'])
        if dmin == 0:
            dprime = 0
        else:
            dprime = d / dmin
    elif round(d, 6) == 0:
        dprime = 0

    # r2
    if values['p1'] == 0 or values['q1'] == 0 or values['p2'] == 0 or values['q2'] == 0:
        r2 = 0
    else:
        dsquared = d**2
        r2 = dsquared/(values['p1'] * values['q1'] * values['p2'] * values['q2'])

    return d, dprime, r2
            

def all_ld_calc(vcf_file, regions, outfile):
    with open(outfile, 'w') as f:
        f.write('chrom1 pos1 chrom2 pos2 d dprime r2\n')

        ref, target = regions
        ref_sites = snppuller(vcf_file, chrom = ref)

        for record1 in tqdm(ref_sites):
            target_sites = snppuller(vcf_file, chrom = target)
            for record2 in target_sites:
                if get_freqs(record1, record2):
                    d, dprime, r2 = ld_calc(*get_freqs(record1, record2))
                    line_out = [ref, record1.POS, target, record2.POS,
                                d, dprime, r2]
                    f.write(' '.join([str(i) for i in line_out]) + '\n')
                

def main():
    vcf_file, regions, outfile = args()
    all_ld_calc(vcf_file, regions, outfile)

if __name__ == '__main__':
    main()

        

