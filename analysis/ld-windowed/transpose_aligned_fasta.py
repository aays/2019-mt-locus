'''
transpose_aligned_fasta.py -
will take in an input fasta and put it into 'long' format

also creates a binary column to ID SNPs that are suitable
for LD calculations

positions are 1-based
'''

import argparse
from Bio import SeqIO
from tqdm import tqdm


def args():
    parser = argparse.ArgumentParser(description = 'Transpose an input fasta',
                                     usage = 'python3.5 transpose_aligned_fasta.py [options]')

    parser.add_argument('-f', '--fasta', required = True,
                        type = str, help = 'Input aligned FASTA')
    parser.add_argument('-o', '--outfile', required = True,
                        type = str, help = 'Name of file to write to')
    parser.add_argument('-s', '--offset', required = False,
                        type = str, help = 'Value to offset position column by [optional]')
    parser.add_argument('-x', '--snps_only', required = False,
                        action = 'store_true', help = 'Will output usable SNPs only [optional]')

    args = parser.parse_args()

    return args.fasta, args.outfile, int(args.offset), args.snps_only


def snp_check(base_list):
    '''(list) -> str
    checks whether a given position is a usable SNP.
    to return 1, SNPs have to be
        1. variant sites (obviously)
        2. diallelic
        3. non-singletons
    '''
    # remove gaps/Ns from consideration
    base_list = [base for base in base_list if base not in ['-', 'N']]

    is_variant = False
    diallelic = False
    non_singleton = False
    base_set = set(base_list)

    # variant check
    if len(base_set) > 1:
        is_variant = True

    # diallelic
    if len(base_set) == 2:
        diallelic = True

    # non singleton
    base_counts = [base_list.count(base) for base in list(base_set)]
    if 1 not in base_counts:
        non_singleton = True

    if False not in [is_variant, diallelic, non_singleton]:
        return '1'
    else:
        return '0'


def transpose_fastas(infile, outfile, offset = None, snps_only = False):
    ''' (str, str, int, bool) -> None
    takes in aligned fasta and transposes to file
    offsets by position if specified
    will also specify whether positions contain usable SNPs
    '''
    with open(outfile, 'w') as f:
        seqs = [s for s in SeqIO.parse(infile, 'fasta')]
        strains = [s.id for s in seqs]
        seq_len = len(seqs[0].seq)
        assert len(set([len(s.seq) for s in seqs])) == 1 # all same length

        # header
        f.write('position ' + ' '.join(strains) + ' is_snp\n') 
        
        for i in tqdm(range(seq_len)):
            if offset:
                position = i + offset + 1
            elif not offset:
                position = i
            bases = [seq[i] for seq in seqs]
            bases_out = ' '.join([seq[i] for seq in seqs])
            is_snp_out = snp_check(bases)

            if snps_only:
                if is_snp_out == '0':
                    continue
                elif is_snp_out == '1':
                    f.write(str(position) + ' ' + bases_out + ' ' + is_snp_out + '\n')
            elif not snps_only:
                f.write(str(position) + ' ' + bases_out + ' ' + is_snp_out + '\n')



def main():
    fasta, outfile, offset, snps_only = args()
    transpose_fastas(fasta, outfile, offset, snps_only)

if __name__ == '__main__':
    main()

        

