'''
append_c_domain.py - adds C domain to main alignment
'''

import argparse
import re
from Bio import SeqIO
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(description = 'Adds C domain at end of main alignment',
                                     usage = 'python3.5 append_c_domain.py [options]')

    parser.add_argument('-m', '--mt_aligned', required = True,
                        type = str, help = 'FASTA containing aligned mt')
    parser.add_argument('-c', '--c_domain', required = True,
                        type = str, help = 'FASTA containing C domain')
    parser.add_argument('-o', '--outfile', required = True,
                        type = str, help = 'File to write to.')

    args = parser.parse_args()

    return args.mt_aligned, args.c_domain, args.outfile


def add_c_domain(mt_aligned, c_domain, outfile):
    strains = [s.id for s in SeqIO.parse(mt_aligned, 'fasta')]
    seqs_out = dict.fromkeys(strains, '')

    # add mt aligned sequences
    for record in SeqIO.parse(mt_aligned, 'fasta'):
        seqs_out[record.id] += str(record.seq)

    # add C domain sequences
    for record in tqdm(SeqIO.parse(c_domain, 'fasta')):
        id_corrected = re.search('CC[0-9]{4}', record.id).group(0)
        seqs_out[id_corrected] += str(record.seq)

    with open(outfile, 'w') as f:
        for strain in seqs_out:
            f.write('>' + strain + '\n')
            f.write(seqs_out[strain] + '\n')
    
def main():
    mt_aligned, c_domain, outfile = args()
    add_c_domain(mt_aligned, c_domain, outfile)

if __name__ == '__main__':
    main()

        

