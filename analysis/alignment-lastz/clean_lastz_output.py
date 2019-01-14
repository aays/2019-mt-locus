'''
clean_lastz_output-2.py - part 2 of the lastz output filtering scripts

if a target region has multiple overlapping hits, this script
will resolve them by picking the highest scoring one

takes in a threshold parameter - if threshold is
0.7, then if alignment B is overlapping with alignment A
such that the overlap is >=70% of len(alignment B), 
then alignment B is also filtered out
'''

import argparse
import csv
from copy import deepcopy
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(description = 'clean_lastz_output-2.py',
                                     usage = 'python3.5 clean_lastz_output-2.py [options]')

    parser.add_argument('-f', '--filename', required = True,
                        type = str, help = 'LASTZ infile (format=general)')
    parser.add_argument('-t', '--threshold', required = True,
                        type = float, help = 'Overlap threshold [0-1].')
    parser.add_argument('-o', '--outfile', required = True,
                        type = str, help = 'File to write to.')

    args = parser.parse_args()

    return args.filename, args.threshold, args.outfile

def in_interval(line1, line2, threshold):
    '''(list, list) -> bool
    takes in two lists of size 3,
    where each list contains:
    1. alignment score
    2. start
    3. end
    and returns whether list B:
    1. overlaps with list A over a given threshold
        - ie if >{threshold}% of region B is in A
    2. has a lower score than list A

    helper function for create_lookup
    '''
    score1, start1, end1 = line1
    score2, start2, end2 = line2
    overlap = len(
        set(range(start1, end1)).intersection(
        set(range(start2, end2)))
    ) / (end2 - start2)

    if start1 < end1 < start2 < end2:
        return False

    elif overlap < threshold:
        return False # keep both regardless

    # aln2 is subset of aln1 and aln1 is better
    elif start1 <= start2 and end1 >= end2 \
    and score1 > score2 and overlap >= threshold:
        return 'left'

    # aln2 overlaps right boundary of aln1
    elif start1 <= start2 < end1 <= end2 and overlap >= threshold:
        if score1 > score2:
            return 'left'
        elif score1 < score2:
            return 'right'
        
def create_lookup(filename, threshold):
    with open(filename, 'r') as f:
        reader = csv.DictReader(f, delimiter = '\t')
        lines = [[int(l['#score']), int(l['zstart1']), int(l['end1'])]
                 for l in reader]
    # reduce lines
    filtered_lines = deepcopy(lines)
    for line1 in lines:
        for line2 in lines:
            overlap_status = in_interval(line1, line2, threshold)
            if overlap_status == 'left':
                if line2 in filtered_lines:
                    filtered_lines.remove(line2)
            elif overlap_status == 'right':
                try:
                    if line1 in filtered_lines:
                        filtered_lines.remove(line1)
                except:
                    print(line1, line2)
            elif overlap_status == False:
                continue
    return filtered_lines

def get_fieldnames(filename):
    with open(filename, 'r') as f:
        fieldnames = f.readline().rstrip('\n').split('\t')
    return fieldnames

def filter_output(filename, filtered_lines, outfile):
    with open(outfile, 'w') as f_out:
        fieldnames = get_fieldnames(filename)
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        with open(filename, 'r') as f_in:
            reader = csv.DictReader(f_in, delimiter='\t')
            for record in tqdm(reader):
                test_parameters = [int(num) for num in
                [record['#score'], record['zstart1'], record['end1']]]
                if test_parameters in filtered_lines:
                    writer.writerow(record)
            

def main():
    filename, threshold, outfile = args()
    filtered_lines = create_lookup(filename, threshold)
    filter_output(filename, filtered_lines, outfile)

if __name__ == '__main__':
    main()

        
