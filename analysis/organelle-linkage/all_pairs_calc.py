'''
all_pairs_calc.py - calculate LD stats between all regions in a VCF
'''

import argparse
import vcf
from tqdm import tqdm
from ld_calc import get_freqs, ld_calc

def args():
    parser = argparse.ArgumentParser(description = 'all_pairs_calc.py',
                                     usage = 'python3.5 all_pairs_calc.py [options]')

    parser.add_argument('-v', '--vcf_file', required = True,
                        type = str, help = 'VCF file to calculate LD from.')
    parser.add_argument('-o', '--outfile', required = True,
                        type = str, help = 'File to write to.')

    args = parser.parse_args()

    return args.vcf_file, args.outfile


def snppuller(vcf_file):
    '''(str) -> iterable
    creates an iterable that lazily loads in usable SNPs.
    usable SNPs must be diallelic and nonsingletons.

    this version of snppuller differs from the ld_calc
    one in that it doesn't take in chrom, start, or end values
    '''
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

    for record in vcfin:
        if usable_record(record):
            yield record
        else:
            continue

def straingetter(record1, record2, GQ_threshold = 30, DP_threshold = 20):
    '''(vcf.Record, vcf.Record, int, int)
    returns a list of strains with calls at both sites.
    also filters for GQ and DP thresholds (defaults - GQ = 30, DP = 20)
    '''
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

    return d, dprime, r2

def all_pairs_ld_calc(vcf_file, outfile):
    '''(str, list, str) -> None
    iterates through pairwise combinations of SNPs in input
    file and writes LD values to outfile.
    '''
    with open(outfile, 'w') as f:
        f.write('chrom1 pos1 chrom2 pos2 d dprime r2\n')

        ref_sites = snppuller(vcf_file)

        for record1 in tqdm(ref_sites):
            if len(record1.ALT) > 1:
                continue
            target_sites = snppuller(vcf_file)
            for record2 in target_sites:
                if len(record2.ALT) > 1:
                    continue
                elif get_freqs(record1, record2):
                    d, dprime, r2 = ld_calc(*get_freqs(record1, record2))
                    line_out = [record1.CHROM, record1.POS, record2.CHROM, 
                                record2.POS, d, dprime, r2]
                    f.write(' '.join([str(i) for i in line_out]) + '\n')
                

def main():
    vcf_file, outfile = args()
    all_pairs_ld_calc(vcf_file, outfile)

if __name__ == '__main__':
    main()

        

