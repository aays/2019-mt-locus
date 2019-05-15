'''
zns_calc_tabix.py - calculate Kelly's ZnS in non-overlapping windows

this script expects the r2 input file to be tabix indexed, and will run
faster than zns_calc.py

for the file to be tabix indexed, it will require a chromosome column -
this script can do that prior to zns calculation if needed (--tabix_for_me)
'''

import argparse
import csv
import tabix
import subprocess
import re
import os
import time
import sys
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(description = 'Calculate ZnS in non-overlapping windows, feat. tabix',
                                     usage = 'python3.5 zns_calc_tabix.py [options]')

    parser.add_argument('-f', '--filename', required = True,
                        type = str, help = 'Input LD file from r2_calc.py')
    parser.add_argument('-w', '--windowsize', required = True,
                        type = int, help = 'Windowsize')
    parser.add_argument('-r', '--calc_range', required = True,
                        type = str, help = 'Start-end, in the tabix-like format "x-y"')
    parser.add_argument('-o', '--outfile', required = True,
                        type = str, help = 'File to write to')
    parser.add_argument('-t', '--tabix_for_me', required = False,
                        action = 'store_true', help = 'This option will tabix-index the input file for you')

    args = parser.parse_args()

    return args.filename, args.windowsize, args.calc_range, \
           args.outfile, args.tabix_for_me


def prep_file(filename):
    ''' (str) -> None
    if infile is not tabixed, will:
    1) add a chromosome column (chrs inferred from filename)
    2) bgzip the file with subprocess
    3) tabix the file

    both the bgzipped file and the tabix index will be stored
    in the same dir as the original infile.
    '''
    chrom = re.search('chromosome_[0-9]{1,2}', filename).group()
    print('Preparing input file for tabix indexing...')
    with open(filename + '_temp', 'w') as f_out:
        f_out.write('\t'.join(['chrom', 'snp1', 'snp2', 'r2']) + '\n')
        with open(filename, 'r') as f_in:
            for line in f_in:
                if line.startswith('snp1'):
                    continue
                line = line.split(' ')
                line.insert(0, chrom)
                f_out.write('\t'.join(line))

    print('bgzipping and tabixing...')
    subprocess.call('bgzip {fname}'.format(fname=filename + '_temp'), shell=True)
    time.sleep(3)
    os.rename(filename + '_temp.gz', filename + '.gz')
    subprocess.call('tabix -p vcf {fname}'.format(fname=filename + '.gz'), shell=True)
    time.sleep(3)
    print('File has been tabix indexed.')


def windowed_zns(filename, windowsize, calc_range, outfile):
    '''(str, int, str) -> None
    calculates zns in windows given input file from r2_calc.
    '''
    try:
        start, end = [int(i) for i in calc_range.split('-')]
        if start == 1:
            start = 0
    except:
        print('Coordinates incorrectly formatted.')
        print('Coordinates provided:', calc_range)
        print('Required format for 1-2 Mb: 1000000-2000000')

    # infer chromosome name from file
    chrom = re.search('chromosome_[0-9]{1,2}', filename).group()

    with open(outfile, 'w') as f_out:
        f_out.write('start end zns site_count\n')

        for window_left in tqdm(range(start, end, windowsize)):
            window_r2 = 0.0
            window_right = window_left + windowsize
            sites = []

            f_tabix = tabix.open(filename)
            for record in f_tabix.query(chrom, window_left, window_right):
                snp1, snp2, r2 = [float(x) for x in record[1:]]
                if window_left <= snp1 < window_right and window_left < snp2 < window_right:
                    window_r2 += r2
                    sites.extend([snp1, snp2])
                    sites = list(set(sites))
                    continue
                elif snp2 > window_right:
                    continue
            sites = list(set(sites)) # once more to be sure
            site_count = len(sites)
            if site_count:
                coefficient = 2 / (site_count * (site_count - 1))
                window_zns = window_r2 * coefficient

                line_out = [str(i) for i in
                            [window_left, window_right, window_zns, site_count]]
                f_out.write(' '.join(line_out) + '\n')
                continue
            elif not site_count:
                line_out = [str(i) for i in [window_left, window_right, 0, 0]]
                f_out.write(' '.join(line_out) + '\n')
                continue
                        
def main():
    filename, windowsize, calc_range, outfile, tabix_for_me = args()
    if not tabix_for_me and not filename.endswith('.gz'):
        print('Your input file does not seem to be bgzipped...')
        print('but you did not ask the script to do the bgzip/tabix operation for you')
        print('This can be done with the --tabix_for_me flag if needed. Exiting...')
        sys.exit(1)
    if tabix_for_me:
        prep_file(filename)
        filename += '.gz'
    windowed_zns(filename, windowsize, calc_range, outfile)

if __name__ == '__main__':
    main()

