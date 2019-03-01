'''
zns_genes.py - calculates zns for directory containing genes
'''

import argparse
import re
import os
import csv
import glob
import ast
from Bio import SeqIO
from tqdm import tqdm
from transpose_aligned_fasta import snp_check
from r2_calc import get_freqs, dcalc, r2calc, usable_pair, ld_calc


def args():
    parser = argparse.ArgumentParser(description = 'calculate zns for directory containing genes',
                                     usage = 'python3.5 zns_genes.py [options]')

    parser.add_argument('-f', '--filename', required = True,
                        type = str, help = 'csv file containing gene-specific stats')
    parser.add_argument('-d', '--directory', required = True,
                        type = str, help = 'Directory containing aligned FASTAs')
    # analysis/cds-popgen/VCF2FASTA/mtLimited -- limited genes
    # analysis/cds-popgen/TranslationAligner/curated -- shared genes
    parser.add_argument('-g', '--gene_type', required = True,
                        type = str, help = 'Type of genes [shared|limited]')
    parser.add_argument('-o', '--outfile', required = True,
                        type = str, help = 'File to write to.')

    args = parser.parse_args()

    return args.filename, args.directory, args.gene_type, args.outfile

def transpose_exons(infile, gene_coords, exon_coords, outfile):
    ''' (str, str, str, str) -> None
    takes in aligned fasta + exon coords and transposes to file
    will also specify whether positions contain usable SNPs
    modified from transpose_aligned_fasta.py

    exon reduction: https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-list-of-lists
    '''
    with open(outfile, 'w') as f:
        seqs = [s for s in SeqIO.parse(infile, 'fasta')]
        strains = [re.search('CC[0-9]{4}', s.id).group(0) for s in seqs]
        seq_len = len(seqs[0].seq)
        assert len(set([len(s.seq) for s in seqs])) == 1 # all same length

        # header
        f.write('position ' + ' '.join(strains) + ' is_snp\n') 

        chrom, coords = gene_coords.split(':')
        start, end = [int(num) for num in coords.split('-')]
        exons = ast.literal_eval(exon_coords)
        exons_all = [list(range(start, end)) for start, end in exons]
        exons_all = [pos for sublist in exons_all for pos in sublist]
        
        for i in tqdm(range(seq_len)):
            position = i + start
            if position not in exons_all:
                continue
            else:
                bases = [seq[i] for seq in seqs]
                bases_out = ' '.join([seq[i] for seq in seqs])
                is_snp_out = snp_check(bases)

                f.write(str(position) + ' ' + bases_out + ' ' + is_snp_out + '\n')

def create_transposed_reference(gene_name, gene_coords, exon_coords, directory):
    ''' (str, str, str, str) -> None
    path-aware wrapper for transpose_exons above that creates a temp dir
    and writes transposed fastas to it
    '''
    if not os.path.exists('temp_transposed'):
        os.mkdir('temp_transposed')
    if not directory.endswith('/'):
        directory += '/'

    # get reference fasta
    pattern = directory + gene_name + '*' + 'fasta'
    globbed = glob.glob(pattern)
    assert len(globbed) == 1
    fname = globbed[0]

    # create a temp transposed file w/ exons only
    outfile = 'temp_transposed/' + gene_name + '_transposed.txt'
    transpose_exons(fname, gene_coords, exon_coords, outfile)
    

def get_gene_zns(gene_name, gene_type):
    ''' (str, str) -> float
    uses the ld_calc function from the r2_calc script to generate
    temporary r2 files using the transposed fastas created above.
    once done, will calculate zns across CDSs of that gene.
    '''
    fname = 'temp_transposed/' + gene_name + '_transposed.txt'

    # get length of sequence for 'windowsize' parameter
    with open(fname, 'r') as f:
        seq_len = len([line for line in f.readlines()]) - 1 # remove header

    # write temp r2 file
    r2_file = 'temp_transposed/' + gene_name + '_r2.txt'

    # maybe 'non_gametolog' wasn't the best choice of words...
    if gene_type == 'shared':
        non_gametolog = False
    elif gene_type == 'limited':
        non_gametolog = True
    ld_calc(fname, r2_file, windowsize=seq_len, non_gametolog=non_gametolog)

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


def all_genic_zns(filename, directory, gene_type, outfile):
    ''' (str, str, str, str) -> None
    'master' function that calculates zns for each gene above
    and outputs a new csv containing a zns_all column.
    '''
    with open(outfile, 'w') as f_out:
        with open(filename, 'r') as f_in:
            records = [record for record in csv.DictReader(f_in)]
            fieldnames = sorted(list(records[0].keys()))
            fieldnames.append('zns_all') # add zns header
            writer = csv.DictWriter(f_out, fieldnames=fieldnames)
            writer.writeheader()

            for record in records:
                if gene_type == 'limited':
                    dict_out = record
                    if record['mtPlus_coords'] and not record['mtMinus_coords']: # plus gene
                        create_transposed_reference(record['gene'], record['mtPlus_coords'], 
                                                    record['CDSs'], directory)
                    elif record['mtMinus_coords'] and not record['mtPlus_coords']: # minus gene
                        create_transposed_reference(record['gene'], record['mtMinus_coords'], 
                                                    record['CDSs'], directory)
                elif gene_type == 'shared':
                    # using plus coordinates as reference here
                    dict_out = record
                    create_transposed_reference(record['gene'], record['mtPlus_coords'],
                                                record['mtPlus_CDSs'], directory)
                dict_out['zns_all'] = get_gene_zns(record['gene'], gene_type)
                writer.writerow(dict_out)


def main():
    filename, directory, gene_type, outfile = args()
    all_genic_zns(filename, directory, gene_type, outfile)

if __name__ == '__main__':
    main()

        

