'''
mask_paralogs.py - mask user-defined paralogous regions in FASTA file

for use with nongametologous FASTAs that have 'false variants'
introduced due to paralogy
'''

import argparse
import re
from Bio import SeqIO
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(description = '',
                                     usage = 'python3.5 script.py [options]')

    parser.add_argument('-f', '--filename', required = True,
                        type = str, help = 'FASTA to mask')
    parser.add_argument('-m', '--mask_intervals', required = True,
                        type = str, help = 'Text file containing regions to mask (1-based, SAM format)')
    parser.add_argument('-s', '--offset', required = False,
                        type = int, help = 'Offset positions in intervals file?')
    parser.add_argument('-o', '--outfile', required = True,
                        type = str, help = 'File to write to.')

    args = parser.parse_args()

    return args.filename, args.mask_intervals, args.offset, args.outfile


def parse_mask_intervals(mask_intervals, offset):
    '''(str, int) -> list
    the offset setting is useful if the FASTA doesn't start
    at 'position 1', but the intervals do
    '''
    with open(mask_intervals, 'r') as f:
        coordinates = [line.rstrip('\n').split(':')[1] for line in f.readlines()]
    coordinates = [[int(i) for i in line.split('-')] for line in coordinates]
    regions = []
    for start, end in coordinates:
        regions.extend(list(range(start, end + 1)))
    if offset:
        regions = [num - offset for num in regions]
    return sorted(regions)

def get_fasta_dict(filename):
    strains = [s.id for s in SeqIO.parse(filename, 'fasta')]
    seqs = dict.fromkeys(strains, '')
    for record in SeqIO.parse(filename, 'fasta'):
        seqs[record.id] += str(record.seq)
    return seqs

def get_sequence_length(filename):
    for record in SeqIO.parse(filename, 'fasta'):
        seq_len = len(record.seq)
        break
    return seq_len

def mask_fasta(filename, mask_intervals, offset, outfile):
    seqs = get_fasta_dict(filename)
    seq_len = get_sequence_length(filename)
    seqs_out = dict.fromkeys(list(seqs.keys()), '')
    regions_to_mask = parse_mask_intervals(mask_intervals, offset)
    print(len(seqs), 'strains in file.')
    print('Masking regions...')
    for strain in seqs:
        for position in tqdm(range(seq_len)):
            if position not in regions_to_mask:
                seqs_out[strain] += seqs[strain][position]
                continue
            elif position in regions_to_mask:
                seqs_out[strain] += 'N'
                continue
    with open(outfile, 'w') as f:
        print('Writing masked regions to file...')
        for strain in seqs_out:
            strain_id_cleaned = re.search('CC[0-9]{4}', strain).group(0)
            f.write('>' + strain_id_cleaned + '\n')
            f.write(seqs_out[strain] + '\n')
        print('Done.')

def main():
    filename, mask_intervals, offset, outfile = args()
    mask_fasta(filename, mask_intervals, offset, outfile)

if __name__ == '__main__':
    main()

        

