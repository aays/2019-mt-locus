'''
gc_content.py - calculate GC content at 4D sites of non gametolog FASTA files
'''

import argparse
import ant
from tqdm import tqdm
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
from collections import OrderedDict

def args():
    parser = argparse.ArgumentParser(description = 'GC content at 4D sites',
                                     usage = 'python3.5 gc_content.py [options]')

    parser.add_argument('-f', '--filename', required = True,
                        type = str, help = 'Gametolog-masked FASTA file')
    parser.add_argument('-a', '--annotation', required = True,
                        type = str, help = 'Annotation table')
    parser.add_argument('-w', '--windowsize', required = True,
                        type = int, help = 'Windowsize')
    parser.add_argument('-r', '--region', required = True,
                        type = str, help = 'SAM format coordinates of FASTA [1-based]')
    parser.add_argument('-o', '--outfile', required = True,
                        type = str, help = 'File to write to')

    args = parser.parse_args()

    return args.filename, args.annotation, args.windowsize, \
           args.region, args.outfile

def get_consensus(filename):
    ''' (str) -> str
    reads in input aligned fasta and generates a consensus
    consensus uses AlignIO's default threshold of 70%
    '''
    # https://www.biostars.org/p/284637/
    alignment = AlignIO.read(filename, 'fasta')
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.dumb_consensus(threshold=0.7, ambiguous='N')
    return str(consensus)

def parse_region(region):
    ''' (str) -> str, int, int
    parses region specification - used for annotation lookup
    '''
    chromosome = region.split(':')[0]
    start, end = [int(num) for num in region.split(':')[1].split('-')]
    return chromosome, start, end

def gc_content_calc(consensus, annotation, windowsize, region, outfile):
    ''' (str, str, int, str, str) -> None
    uses the consensus sequence generated above to calculate GC content
    at selectively unconstrained sites (intronic, intergenic, 4D).
    writes count of GC nucleotides + Ns + all non-N sites to specified outfile.
    '''
    chromosome, start, end = parse_region(region)
    with open(outfile, 'w') as f:
        f.write('start end GC N total_sites\n')
        seq_index = 0
        for window in tqdm(range(start, end, windowsize)):
            counter = OrderedDict.fromkeys(['GC_count', 'N_count', 'total_sites'], 0)
            window_start = window
            if window + windowsize > end:
                window_end = end
            else:
                window_end = window + windowsize
            p = ant.Reader(annotation)
            for record in p.fetch(chromosome, window_start, window_end):
                seq_index = record.pos - start # fix offset for string indexing
                if record.is_intronic or record.is_intergenic or record.is_fold4:
                    if consensus[seq_index] in ['G', 'C']:
                        counter['GC_count'] += 1
                        counter['total_sites'] += 1
                    elif consensus[seq_index] in ['A', 'T']: # Ns not counted
                        counter['total_sites'] += 1
                    elif consensus[seq_index] == 'N':
                        counter['N_count'] += 1
            window_out = ' '.join([str(num) for num in [window_start, window_end]])
            line_out_counts = ' '.join([str(i) for i in list(counter.values())])
            f.write(window_out + ' ' + line_out_counts + '\n')


def main():
    filename, annotation, windowsize, region, outfile = args()
    print('Obtaining consensus...')
    consensus = get_consensus(filename)
    print('Done.')
    print('Calculating GC content...')
    gc_content_calc(consensus, annotation, windowsize, region, outfile)

if __name__ == '__main__':
    main()

        
