'''
zns_calc.py - calculate Kelly's ZnS in non-overlapping windows
'''

import argparse
import csv
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(description = 'Calculate ZnS in non-overlapping windows',
                                     usage = 'python3.5 zns_calc.py [options]')

    parser.add_argument('-f', '--filename', required = True,
                        type = str, help = 'Input LD file from r2_calc.py')
    parser.add_argument('-w', '--windowsize', required = True,
                        type = int, help = 'Windowsize')
    parser.add_argument('-o', '--outfile', required = True,
                        type = str, help = 'File to write to')

    args = parser.parse_args()

    return args.filename, args.windowsize, args.outfile

def get_file_stats(filename, windowsize):
    '''(str, int) -> int, int, int
    iterates through infile to determine windows
    for windowed zns calculation
    '''
    line_count = 0
    with open(filename, 'r') as f:
        for line in f:
            if line_count == 1:
                start = int(line.split(' ')[0])
            line_count += 1

    with open(filename, 'r') as f:
        counter = 0
        for line in f:
            counter += 1
            if counter == line_count:
                end = int(line.split(' ')[1])

    start = round(start / windowsize) * windowsize
    end = (round(end / windowsize) * windowsize) + windowsize
    return start, end, line_count

def windowed_zns(filename, windowsize, outfile):
    '''(str, int, str) -> None
    calculates zns in windows given input file from r2_calc.
    '''
    with open(outfile, 'w') as f_out:
        start, end, line_count = get_file_stats(filename, windowsize)
        f_out.write('start end zns site_count\n')

        for window_left in tqdm(range(start, end, windowsize)):
            window_r2 = 0.0
            window_right = window_left + windowsize
            sites = []
            with open(filename, 'r') as f_in:
                reader = csv.DictReader(f_in, delimiter = ' ')
                for record in reader:
                    if int(record['snp1']) < window_left:
                        continue

                    elif window_left <= int(record['snp2']) < window_right:
                        if window_left <= int(record['snp1']) < window_right:
                            window_r2 += float(record['r2'])
                            sites.extend([record['snp1'], record['snp2']])
                            sites = list(set(sites)) # reduce redundancy
                            continue

                    elif int(record['snp2']) > window_right:
                        if int(record['snp1']) < window_right:
                            continue
                        elif int(record['snp1']) > window_right:
                            sites = list(set(sites)) # once more to be sure
                            site_count = len(sites)
                            if site_count:
                                coefficient = 2 / (site_count * (site_count - 1))
                                window_zns = window_r2 * coefficient

                                line_out = [str(i) for i in 
                                            [window_left, window_right, window_zns, site_count]]
                                f_out.write(' '.join(line_out) + '\n')
                                break
                            elif not site_count:
                                line_out = [str(i) for i in [window_left, window_right, 0, 0]]
                                f_out.write(' '.join(line_out) + '\n')
                                break

                continue
                        
def main():
    filename, windowsize, outfile = args()
    windowed_zns(filename, windowsize, outfile)

if __name__ == '__main__':
    main()

