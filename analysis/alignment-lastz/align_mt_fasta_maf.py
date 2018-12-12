'''
align_mt_fasta_maf.py - create aligned fastas from maf file
'''

from tqdm import tqdm
import argparse
from Bio.Seq import Seq
from Bio import SeqIO

def args():
    parser = argparse.ArgumentParser(description = 'Create fasta file containing mt+ and mt- sequences',
                                     usage = 'align_mt_fasta_maf.py [options]')

    parser.add_argument('-p', '--plus', required = True,
                        type = str, help = 'FASTA file containing mt+ sequences.')
    parser.add_argument('-m', '--minus', required = True,
                        type = str, help = 'FASTA file containing mt- sequences.')
    parser.add_argument('-a', '--alignment', required = True,
                        type = str, help = 'LASTZ alignment output (--format=maf).')
    parser.add_argument('-b', '--bed', required = True,
                        type = str, help = 'LASTZ alignment output (--format=bed).')
    parser.add_argument('-o', '--outdir', required = True,
                        type = str, help = 'Directory to write to.')

    args = parser.parse_args()

    return [args.plus, args.minus, args.alignment, args.outdir]

class aln_bed(object):
    '''
    Quick class to parse lastz alignment output (bed).
    Expects --format=general formatted output from command-line lastz
    '''
    def __init__(self, score, name1, strand1, size1, zstart1, 
        end1, name2, strand2, size2, zstart2, end2, identity, 
        idPct, coverage, covPct):
        self.score = score
        self.name1 = name1
        self.strand1 = strand1
        self.size1 = size1
        self.zstart1 = int(zstart1) # origin-zero
        self.end1 = int(end1)
        self.name2 = name2
        self.strand2 = strand2
        self.size2 = size2
        self.zstart2 = int(zstart2) # origin-zero, orientation dependent
        self.end2 = int(end2)
        self.identity = [int(num) for num in str(identity).split('/')]
        self.idPct = float(idPct[:len(idPct) - 1])
        self.coverage = [int(num) for num in str(coverage).split('/')]
        self.covPct = covPct[:len(covPct) - 1]

class aln(object):
    '''
    Quick class to parse maf lastz alignment output.
    Expects --format=maf formatted output from command-line lastz
    '''
    def __init__(self, score, name1, start1, match_size1, strand1,
        size1, seq1, name2, start2, match_size2, strand2, size2, seq2):
        self.score = score
        self.name1 = name1
        self.start1 = int(start1) # origin-zero
        self.match_size1 = int(match_size1)
        self.strand1 = strand1
        self.size1 = size1
        self.seq1 = SeqIO.SeqRecord(Seq(seq1.rstrip('\n')))
        self.name2 = name2
        self.start2 = int(start2) # origin-zero, orientation dependent
        self.match_size2 = int(match_size2)
        self.strand2 = strand2
        self.size2 = size2
        self.seq2 = SeqIO.SeqRecord(Seq(seq2.rstrip('\n')))

class aln_file(object):
    def __init__(self, filename):
        with open(filename, 'r') as f:
            d = [line for line in f.readlines()
                 if not line.startswith('#')]

        alignments = []
        chunks = [d[i:i+4]
                  for i in range(0, len(d), 4)]

        for curr_aln in chunks:
            score = curr_aln[0].split('=')[1].rstrip('\n')
            for i in [1,2]: # sequences
                curr_aln[i] = [element for element 
                               in curr_aln[i].split(' ')
                               if element] # not empty

            name1, start1, match_size1, strand1, \
                size1, seq1 = curr_aln[1][1:]
            name2, start2, match_size2, strand2, \
                size2, seq2 = curr_aln[2][1:]

            aln_obj = aln(score, name1, start1, match_size1,
                strand1, size1, seq1, name2, start2,
                match_size2, strand2, size2, seq2)

            alignments.append(aln_obj)

        self._alignments = alignments

    def alignments(self):
        return self._alignments

                
def parse_aln(filename):
    '''
    (str) -> list
    parses LASTZ alignment file to return correct starts
    duplicates should have been removed using R script
    expects --format=general
    '''
    with open(filename) as f:
        aln_bed_file = [aln(*line.split('\t')) 
                        for line in f.readlines() 
                        if not line.startswith('score')]
        # choosing best lastz hit
        correct_starts = dict.fromkeys([a.zstart1
                                      for a in aln_bed_file])
        correct_starts = {}
        for a in aln_bed_file:
            correct_starts[a.zstart1] = [a.score, a.zstart2]

    return correct_starts


def get_strain_refs(filename, reverse = False):
    seqs = {}
    if not reverse:
        for record in SeqIO.parse(filename, 'fasta'):
            seqs[record.id] = record.seq
    elif reverse:
        for record in SeqIO.parse(filename, 'fasta'):
            seqs[record.id] = record.reverse_complement().seq
    return seqs

        

def write_fastas(aln_file, correct_starts, plus_dict, minus_dict, minus_rev_dict, outdir):
    for a in aln_file.alignments():
        score, start1, start2 = a.score, a.start1, a.start2
        # check to see this is the correct alignment
        if correct_starts[start1] != [score, start2]:
            continue
        elif correct_starts[start1] == [score, start2]:
            plus_seqs_out = dict.fromkeys(list(plus_seqs.keys()), '')
            minus_seqs_out = dict.fromkeys(list(minus_seqs.keys()), '')
            chr6_start = 298298 + start1
            chr6_end = 298298 + start1 + a.size1

            # plus sequence
            for i, base in enumerate(a.seq1):
                for strain in plus_seqs_out:
                    if base == '-':
                        plus_seqs_out[strain] += '-'
                    else:
                        plus_seqs_out[strain] += plus_dict[strain][i]

            # minus sequence
            for i, base in enumerate(a.seq2):
                if a.strand2 == '+':
                    for strain in minus_seqs_out:
                        if base == '-':
                            minus_seqs_out[strain] += '-'
                        else:
                            minus_seqs_out[strain] += minus_dict[strain][i]
                elif a.strand2 == '-':
                    for strain in minus_seqs_out:
                        if base == '-':
                            minus_seqs_out[strain] += '-'
                        else:
                            minus_seqs_out += minus_rev_dict[strain][i]

            outname = outdir + 'chromosome_6_' + str(chr6_start) + str(chr6_end)
            seq_id_plus = 'chromosome_6:' + str(chr6_start) + '-' + str(chr6_end)
            seq_id_minus = 'mtMinus:' + str(a.start2) + '-' + \
                str(a.start2 + a.size2) + '|orientation=' + str(a.strand2)
            with open(outname, 'w') as f:
                for strain in plus_seqs_out.keys():
                    f.write('>' + seq_id_plus + '\n')
                    f.write(plus_seqs_out[strain] + '\n')
                for strain in minus_seqs_out.keys():
                    f.write('>' + seq_id_minus + '\n')
                    f.write(minus_seqs_out[strain] + '\n')

def main():
    plus, minus, alignment, bed, outdir = args()

    # parse aln files
    print('Parsing alignment files...')
    aln_maf = aln_file(alignment)
    correct_starts = parse_aln(bed)

    # get refs
    print('Parsing FASTA files...')
    plus_refs = get_strain_refs(plus)
    minus_refs = get_strain_refs(minus)
    minus_rev_refs = get_strain_refs(minus, reverse = True)

    # write fastas
    print('Writing aligned FASTA files...')
    write_fastas(aln_maf, correct_starts, 
        plus_refs, minus_refs, minus_rev_refs, outdir)

    print('Done.')

if __name__ == '__main__':
    main()




                


        












