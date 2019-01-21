'''
zns_genes.py - calculates zns for directory containing genes
'''

import argparse
import re
import os
import csv
import glob
from Bio import SeqIO
from tqdm import tqdm
from transpose_aligned_fasta import snp_check, transpose_fastas
from r2_calc import get_freqs, dcalc, r2calc, usable_pair, ld_calc


def args():
    parser = argparse.ArgumentParser(description = 'calculate zns for directory containing genes',
                                     usage = 'python3.5 zns_genes.py [options]')

    parser.add_argument('-f', '--filename', required = True,
                        type = str, help = 'csv file containing gene-specific stats')
    parser.add_argument('-d', '--directory', required = True,
                        type = str, help = 'Directory containing aligned FASTAs')
    parser.add_argument('-o', '--outfile', required = True,
                        type = str, help = 'File to write to.')

    args = parser.parse_args()

    return args.filename, args.directory, args.outfile

def create_transposed_reference(gene_name, directory):
    if not os.path.exists('temp_transposed'):
        os.mkdir('temp_transposed')
    if not directory.endswith('/'):
        directory += '/'

    # get reference fasta
    pattern = directory + gene_name + '*' + 'fasta'
    globbed = glob.glob(pattern)
    assert len(globbed) == 1
    fname = globbed[0]

    # use transpose_fastas from transpose_aligned_fasta.py
    # this creates a temp transposed file
    outfile = 'temp_transposed/' + gene_name + '_transposed.txt'
    transpose_fastas(fname, outfile, offset=None, snps_only=False)
    

def get_gene_zns(gene_name):
    fname = 'temp_transposed/' + gene_name + '_transposed.txt'

    # get length of sequence for 'windowsize' parameter
    with open(fname, 'r') as f:
        seq_len = len([line for line in f.readlines()]) - 1 # remove header

    # write temp r2 file
    r2_file = 'temp_transposed/' + gene_name + '_r2.txt'
    ld_calc(fname, r2_file, windowsize=seq_len, non_gametolog=False)

    r2_cumulative = 0.0
    sites = []
    with open(r2_file, 'r') as f:
        reader = csv.DictReader(f, delimiter=' ')
        for record in reader:
            r2_cumulative += float(record['r2'])
            sites.extend([record['snp1'], record['snp2']])
            sites = list(set(sites))

    sites = list(set(sites)) # once more to be sure
    site_count = len(sites)
    if site_count:
        coefficient = 2 / (site_count * (site_count - 1)) 
        zns_out = r2_cumulative * coefficient
    elif site_count == 0:
        zns_out = None

    return zns_out


def all_genic_zns(filename, directory, outfile):
    with open(outfile, 'w') as f_out:
        with open(filename, 'r') as f_in:
            records = [record for record in csv.DictReader(f_in)]
            fieldnames = sorted(list(records[0].keys()))
            fieldnames.append('zns_all') # add zns header
            writer = csv.DictWriter(f_out, fieldnames=fieldnames)
            writer.writeheader()

            for record in records:
                dict_out = record
                create_transposed_reference(record['gene'], directory)
                dict_out['zns_all'] = get_gene_zns(record['gene'])
                writer.writerow(dict_out)
        


def main():
    filename, directory, outfile = args()
    all_genic_zns(filename, directory, outfile)

if __name__ == '__main__':
    main()

        

