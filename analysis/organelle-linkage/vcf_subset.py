'''
vcf_subset.py - subsample input VCFs

can be used to pull specific regions from a VCF or
subsample x% of records in VCF
'''

import argparse
import vcf
import random
import sys
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(description='Filter/subset VCFs',
                                     usage='python3.5 vcf_subset.py [options]')

    parser.add_argument('-v', '--vcf', required=True,
                        type=str, help='Input VCF')
    parser.add_argument('-c', '--chrom', required=False,
                        type=str, nargs='+', help = 'Chromosomes to keep. [Optional]')
    parser.add_argument('-f', '--filter_fraction', required=False,
                        type=float, help='Fraction of sites to keep. [Optional]')
    parser.add_argument('-o', '--outfile', required=True,
                        type=str, help='File to write to [.vcf]')

    args = parser.parse_args()

    return args.vcf, args.chrom, args.filter_fraction, args.outfile


def subset_vcf(vcf_file, chrom, filter_fraction, outfile):
    '''(str, list, float, str) -> dict (or None)
    subsets VCF based on given list of chroms and/or fraction to
    filter down to. if filtering, will return dict reporting records
    iterated over + records actually kept.
    '''
    if filter_fraction:
        counter = dict.fromkeys(['kept', 'total'], 0)
    with open(outfile, 'w') as f:
        if not filter_fraction and not chrom:
            print('Error - no subsetting defined.')
            print('Why run this script then?')
            sys.exit(1)
        elif filter_fraction and not chrom:
            reader = vcf.Reader(filename = vcf_file, compressed = True)
            writer = vcf.Writer(f, reader)
            for record in tqdm(reader):
                counter['total'] += 1
                if random.random() <= filter_fraction:
                    counter['kept'] += 1
                    writer.write_record(record)
        elif len(chrom) == 1:
            reader = vcf.Reader(filename = vcf_file, compressed = True)
            region = reader.fetch(chrom = chrom[0])
            writer = vcf.Writer(f, region)
            for record in tqdm(region):
                if filter_fraction:
                    counter['total'] += 1
                    if random.random() <= filter_fraction:
                        counter['kept'] += 1
                        writer.write_record(record)
                elif not filter_fraction:
                    writer.write_record(record)
        elif len(chrom) > 1:
            for c in chrom:
                reader = vcf.Reader(filename = vcf_file, compressed = True)
                region = reader.fetch(chrom = c)
                writer = vcf.Writer(f, region)
                print('Begin chromosome', c)
                for record in tqdm(region):
                    if filter_fraction:
                        counter['total'] += 1
                        if random.random() <= filter_fraction:
                            counter['kept'] += 1
                            writer.write_record(record)
                    elif not filter_fraction:
                        writer.write_record(record)
    if filter_fraction:
        return counter

def main():
    vcf, chrom, filter_fraction, outfile = args()
    if filter_fraction:
        print('Filtering selected.')
        print('Filtering to', filter_fraction)
        counter = subset_vcf(vcf, chrom, filter_fraction, outfile)
        print('Subsetting complete.')
        print(counter['kept'], 'records kept of', counter['total'])
        print('File written to', outfile)
    elif not filter_fraction:
        print('Subsetting selected.')
        subset_vcf(vcf, chrom, filter_fraction, outfile)
        print('Complete.')
        print('File written to', outfile)

if __name__ == '__main__':
    main()

        

