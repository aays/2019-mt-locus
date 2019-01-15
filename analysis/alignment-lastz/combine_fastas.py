'''
combine_fastas.py - create aligned fasta using transposed file

can also mask crappy regions (ie in case of paralogous alignment)
'''

import argparse
import csv
from collections import OrderedDict
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(description = 'combine_fastas.py',
                                     usage = 'python combine_fastas.py [options]')

    parser.add_argument('-f', '--file', required = True,
                        type = str, help = 'Filtered transposed FASTA file.')
    parser.add_argument('-m', '--mask_intervals', required = False,
                        type = str, help = 'Text file containing regions to mask (1-based, SAM format) [optional]')
    parser.add_argument('-o', '--outfile', required = True,
                        type = str, help = 'Filename to write to.')

    args = parser.parse_args()

    return args.file, args.mask_intervals, args.outfile

def get_strain_names(infile):
    with open(infile, 'r') as f:
        cols = f.readline()
    strains = cols.rstrip('\n').split(' ')[1:]
    return strains

def parse_mask_intervals(mask_intervals):
    with open(mask_intervals, 'r') as f:
        coordinates = [line.rstrip('\n').split(':')[1] for line in f.readlines()]
    coordinates = [[int(i) for i in line.split('-')] for line in coordinates]
    regions = []
    for start, end in coordinates:
        regions.extend(list(range(start, end + 1)))
    return sorted(regions)

def write_fasta(infile, outfile, strains, mask_intervals):
    seqs = OrderedDict.fromkeys(strains, '')
    with open(infile, 'r') as f:
        reader = csv.DictReader(f, delimiter = ' ')
        counter = 298299 # start of mt locus, origin one
        end = 826738 # mt locus ends at 826737
        if mask_intervals:
            mask_regions = parse_mask_intervals(mask_intervals)
        for line in tqdm(reader):
            assert counter < end
            position = int(line['position'])
            if position == counter and not mask_intervals:
                for strain in strains:
                    seqs[strain] += line[strain]
                counter += 1
                continue
            elif position == counter and mask_intervals:
                if position not in mask_regions:
                    for strain in strains:
                        seqs[strain] += line[strain]
                    counter += 1
                    continue
                elif position in mask_regions:
                    for strain in strains:
                        seqs[strain] += 'N'
                    counter += 1
                    continue
            elif position > counter: # gap
                while counter < position:
                    for strain in strains:
                        seqs[strain] += 'N'
                    counter += 1
                if counter == position:
                    for strain in strains:
                        seqs[strain] += line[strain]
                    counter += 1
                    continue
                elif counter > position:
                    print('what the shit!')
    with open(outfile, 'w') as f:
        for strain in seqs.keys():
            f.write('>' + strain + '\n')
            f.write(seqs[strain] + '\n')
            
def main():
    infile, mask_intervals, outfile = args()
    strains = get_strain_names(infile)
    write_fasta(infile, outfile, strains, mask_intervals)

if __name__ == '__main__':
    main()

        

