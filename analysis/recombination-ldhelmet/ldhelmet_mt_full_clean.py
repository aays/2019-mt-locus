'''
ldhelmet_mt_full_clean.py - 
only update coordinates of LDhelmet file
without filtering out any rows

most functions from make_mt_only.py
'''

import argparse
from Bio import SeqIO
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(description = 'ldhelmet_mt_full_clean.py',
                                     usage = 'ldhelmet_mt_full_clean.py [options]')

    parser.add_argument('-i', '--filename', required = True,
                        type = str, help = 'LDhelmet infile')
    parser.add_argument('-a', '--allele', required = True,
                        type = str, help = 'mt allele [plus/minus]')
    parser.add_argument('-o', '--outfile', required = True,
                        type = str, help = 'File to write to')

    args = parser.parse_args()

    return args.filename, args.allele, args.outfile

def correct_line(line, offset):
    ''' (list) -> str
    returns modified line as a str
    '''
    line = line.split(' ')
    line[0] = int(line[0]) + offset
    line[1] = int(line[1]) + offset
    line = [str(i) for i in line]

    return ' '.join(line)

def filter_ldhelmet(infile, allele, outfile):
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
                    line = correct_line(line, offset)
                    out.write(line)

def main():
    filename, allele, outfile = args()

    filter_ldhelmet(filename, allele, outfile)

if __name__ == '__main__':
    main()

