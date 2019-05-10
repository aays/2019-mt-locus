'''
autosomal_ld.py - run entire zns pipeline on an input chromosome fasta

should be run from project root - contains some hardcoded paths
'''

import argparse
import os
import subprocess
import time
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(description = 'autosomal_ld.py - get autosomal ZnS',
                                     usage = 'python3.5 autosomal_ld.py [options]')

    parser.add_argument('-f', '--filename', required = True,
                        type = str, help = 'FASTA file')
    parser.add_argument('-c', '--chrom', required = True,
                        type = str, help = 'Chromosome')

    args = parser.parse_args()

    return args.filename, args.chrom

# hardcoded lengths
lengths = {'chromosome_1': 8033585,
'chromosome_2': 9223677,
'chromosome_3': 9219486,
'chromosome_4': 4091191,
'chromosome_5': 3500558,
'chromosome_6': 9023763,
'chromosome_7': 6421821,
'chromosome_8': 5033832,
'chromosome_9': 7956127,
'chromosome_10': 6576019,
'chromosome_11': 3826814,
'chromosome_12': 9730733,
'chromosome_13': 5206065,
'chromosome_14': 4157777,
'chromosome_15': 1922860,
'chromosome_16': 7783580,
'chromosome_17': 7188315}

def transpose(chrom):
    transpose_cmd = 'time python3.5 analysis/ld-windowed/transpose_aligned_fasta.py \
    --fasta data/aligned-fastas/autosomal/{chrom}_all.fasta \
    --outfile data/ld-windowed/autosomal/{chrom}_long.txt \
    --offset 0'

    print('Transposing', chrom, 'FASTA')
    subprocess.call(transpose_cmd.format(chrom=chrom), shell=True)
    print(chrom, 'Transposition complete.')
    time.sleep(3)

def r2_calc(chrom):
    r2_cmd = 'time python3.5 analysis/ld-windowed/r2_calc_region.py \
    --filename data/ld-windowed/autosomal/{chrom}_long.txt \
    --windowsize 1000 \
    --region {chrom}:{start}-{end} \
    --outfile data/ld-windowed/{chrom}_temp/{chrom}_{start}-{end}.txt'

    windows = list(range(0, lengths[chrom], 1000000))
    windows[0] = 1
    windows.append(lengths[chrom])

    print('Starting r2 calculations for', chrom)
    for i in range(len(windows) - 1):
        start, end = windows[i], windows[i + 1]
        print('r2 for {0} {1} {2}'.format(chrom, start, end))
        subprocess.call(r2_cmd.format(chrom=chrom, start=start, end=end), shell=True)
        time.sleep(3)
    
def zns_calc(chrom):
    zns_cmd = 'time python3.5 analysis/ld-windowed/zns_calc_tabix.py \
    --filename data/ld-windowed/{chrom}_temp/{chrom}_{start}-{end}.txt \
    --windowsize 1000 \
    --calc_range {start}-{end} \
    --outfile data/ld-windowed/{chrom}_temp/{chrom}_{start}-{end}.zns \
    --tabix_for_me'

    windows = list(range(0, lengths[chrom], 1000000))
    windows[0] = 1
    windows.append(lengths[chrom])

    print('ZnS calculations for', chrom)
    for i in range(len(windows) - 1):
        start, end = windows[i], windows[i + 1]
        print('zns for {0} {1} {2}'.format(chrom, start, end))
        subprocess.call(zns_cmd.format(chrom=chrom, start=start, end=end), shell=True)
        time.sleep(3)

def combine(chrom):
    combine_cmd = 'time Rscript analysis/ld-windowed/combine_zns_autosomal.R \
    --directory data/ld-windowed/{chrom}_temp \
    --chrom {chrom} \
    --outfile data/ld-windowed/autosomal/{chrom}_final.zns'

    print('Combining ZnS files for', chrom)
    subprocess.call(combine_cmd.format(chrom=chrom), shell=True)
    time.sleep(1)

def main():
    filename, chrom = args()
    temp_dir = 'data/ld-windowed/' + chrom + '_temp'
    if not os.path.isdir(temp_dir):
        os.mkdir(temp_dir)
    time.sleep(1)
    transpose(chrom)
    r2_calc(chrom)
    zns_calc(chrom)
    combine(chrom)
    

if __name__ == '__main__':
    main()

        

