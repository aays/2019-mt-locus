'''
ldhelmet_clean.py - correct coordinates in LDhelmet outputs

filename coordinates are 0-based, but LDhelmet coordinates are 1-based

files need to have been created w/ align_mt_fasta_maf.py and
ldhelmet_init.sh - this script relies on consistent naming of infiles
'''

import argparse
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(description = 'ldhelmet_clean.py - correct coordinates',
                                     usage = 'ldhelmet_clean.py [options]')

    parser.add_argument('-f', '--filename', required = True,
                        type = str, help = 'File to correct')
    parser.add_argument('-o', '--outfile', required = True,
                        type = str, help = 'File to write to')

    args = parser.parse_args()

    return args.filename, args.outfile

def correct_line(line, offset):
    ''' (list) -> str
    returns modified line as a str
    '''
    line[0] = int(line[0]) + offset
    line[1] = int(line[1]) + offset
    line = [str(i) for i in line]

    return ' '.join(line)

def correct_coordinates(filename, outfile):
    ''' (str) -> None
    iterates through LDhelmet file and corrects coordinates
    '''
    start = 298298
    
    colnames = ['left_snp', 'right_snp', 'mean', 
                'p0.025', 'p0.500', 'p0.975\n']
    
    with open(outfile, 'w') as f_out:
        # add column header
        f_out.write(' '.join(colnames))

        with open(filename, 'r') as f_in:
            for line in tqdm(f_in):
                if line.startswith(('#', 'version')):
                    continue
                else:
                    line = line.split(' ')
                    corrected_line = correct_line(line, offset = start)
                    f_out.write(corrected_line)


def main():
    filename, outfile = args()
    correct_coordinates(filename, outfile)

if __name__ == '__main__':
    main()

