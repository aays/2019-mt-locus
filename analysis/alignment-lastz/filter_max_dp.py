'''
usage:
python3.5 filter_max_dp.py [allele] [fname] [outname] [depth multiplier]
'''

import argparse
import vcf
import sys
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(description = 'Filter by max DP per strain',
                                     usage = 'filter_max_dp.py [options]')

    parser.add_argument('-f', '--filename', required = True,
                        type = str, help = 'Input VCF')
    parser.add_argument('-o', '--outfile', required = True,
                        type = str, help = 'Outfile name')
    parser.add_argument('-a', '--allele', required = True,
                        type = str, help = 'mt allele')
    parser.add_argument('-d', '--dp_multiplier', required = True,
                        type = str, help = 'Max DP multiplier to filter by')

    args = parser.parse_args()

    return args.filename, args.outfile, args.allele, int(args.dp_multiplier)

def get_coordinates(allele):
    if allele == 'plus':
        chrom = 'chromosome_6'
        start = 298299
        end = 826737
    elif allele == 'minus':
        chrom = 'mtMinus'
        start = 1
        end = 345555
    return chrom, start, end

def get_strains(filename):
    vcfin = vcf.Reader(filename = filename, compressed = True)
    rec = next(vcfin)
    samples = [call.sample for call in rec]
    return samples


def mean_dp_per_strain(filename, allele):
    print('Calculating mean DP per strain...')
    samples = get_strains(filename)
    chrom, start, end = get_coordinates(allele)
    dp_vals = dict.fromkeys(samples, 0.0)
    record_counts = dict.fromkeys(samples, 0)
    for strain in tqdm(dp_vals.keys()):
        vcfin = vcf.Reader(filename = filename, compressed = True)
        region = vcfin.fetch(chrom, start, end)
        for record in region:
            for i, call in enumerate(record.samples):
                if call['GT'] == '.':
                    continue
                elif call.sample == strain and call['DP']:
                    dp_vals[strain] += call['DP']
                    record_counts[strain] += 1
    mean_dp_vals = dict.fromkeys(samples, 0.0)
    for strain in dp_vals.keys():
        mean_dp_vals[strain] = dp_vals[strain] / record_counts[strain]

    return mean_dp_vals
            

def make_blanked_call(sample, position, keys = 5):
    if keys == 5:
        sample_format = vcf.model.make_calldata_tuple('GT:AD:DP:GQ:PL'.split(':'))
        sampdat = ['.', None, None, None, None]
    elif keys == 3:
        sample_format = vcf.model.make_calldata_tuple('GT:AD:DP'.split(':'))
        sampdat = ['.', None, None]
    name = sample
    site = position
    blanked_call = vcf.model._Call(site, name, sample_format(*sampdat)) # from pyvcf source code
    return blanked_call


def get_len_call(record):
    # some records have 3 fields, not 5
    # needed to check for those
    call = record.samples[0].data
    return len(call)


def write_filtered_vcf(filename, outfile, allele, dp_multiplier):
    print('Filtering VCF...')
    samples = get_strains(filename)
    chrom, start, end = get_coordinates(allele)
    mean_dp_vals = mean_dp_per_strain(filename, allele)
    kept_counter = dict.fromkeys(samples, 0)
    total_counter = dict.fromkeys(samples, 0)
    with open(outfile, 'w') as f:
        vcfin = vcf.Reader(filename = filename, compressed = True)
        writer = vcf.Writer(f, vcfin)
        region = vcfin.fetch(chrom, start, end)
        for record in tqdm(region):
            for i, call in enumerate(record.samples):
                if call['GT'] == '.':
                    continue
                elif not call['DP']: # weird sites w/ calls but no DP
                    continue
                else:
                    call_dp = call['DP']
                    dp_threshold = mean_dp_vals[call.sample] * dp_multiplier
                    total_counter[call.sample] += 1
                    kept_counter[call.sample] += 1
                    if call_dp > dp_threshold:
                        call_len = get_len_call(record)
                        record.samples[i] = make_blanked_call(call.sample, record.POS, call_len)
                        kept_counter[call.sample] -= 1
            writer.write_record(record)
    return kept_counter, total_counter
                    

def main():
    filename, outfile, allele, dp_multiplier = args()
    kept, total = write_filtered_vcf(filename, outfile, allele, dp_multiplier)
    with open('vcf_filter_log.txt', 'w') as f:
        f.write('strain total kept\n')
        for strain in total.keys():
            line_out = [strain, total[strain], kept[strain]]
            line_out = [str(i) for i in line_out]
            f.write(' '.join(line_out) + '\n')
    print('Completed.')
    print('Log written to vcf_filter_log.txt')

if __name__ == '__main__':
    main()


