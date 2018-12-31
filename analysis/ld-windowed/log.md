log for LD analysis across mt locus

## 25/12/2018

- created plink folders in `data` and `analysis` 
- add plink executable to `bin`

overall to do:
1. run plink across mt locus?
2. write script that takes in tabular plink outfile and calculate ZnS in windows

## 28/12/2018

plink actually won't work for this - our aligned mt locus is in fasta format and can't
be converted to a VCF. plink doesn't take in fasta files as input

more detailed:

`transpose_aligned_fasta.py`
- need a script that returns LD (r2) values from fasta files
- to do this - need a long format file once again
    - save all sequences as SeqObjects in a dict
    - file.write header (ie position strain1 strain2 is_snp)
        - the is_snp column will be 1 for non-singleton diallelic SNPs - makes next script easier
    - for position in sequences:
        - outlist = [position, base1, base2, etc]
        - fxn that reads in outlist[1:] and checks whether this is a usable SNP (non-singleton, diallelic)

```bash
time python3.5 analysis/ld-windowed/transpose_aligned_fasta.py \
--fasta data/aligned-fastas/mt_aligned_final.fasta \
--outfile data/ld-windowed/mt_aligned_long.txt \
--offset 298298 # mt plus coordinate
```

## 31/12/2018

`r2_calc.py`
- inputs - infile, outfile, windowsize
- open windowed LD outfile
    - open transposed file w/ csv reader ('ref')
    - for line in transposed ref file:
        - open transposed file *again* w/ csv reader ('target')
        - if ref_snp[is_snp] == 1:
            - for line in transposed target file:
                - if target_snp.position - ref_snp.position is within windowsize:
                    - ld = r2calc(ref_snp, target_snp)
                    - outfile.write(ref_snp.pos, target_snp.pos, ld)
                
```bash
python3.5 analysis/ld-windowed/r2_calc.py \
--filename data/ld-windowed/mt_aligned_long.txt \
--windowsize 10000 \
--outfile data/ld-windowed/mt_aligned_r2.txt
```

to do tomorrow - maybe a modification of the
fasta transposition script that *only* returns
usable SNP sites? could speed up the r2 script, since it's
rather slow...

also, when showing ZnS in a plot - Kelly's paper (1997)
plots each successive point as 'one more SNP away from
a focal site' (Fig 4). but how does Flowers plot ZnS?...

## 31/12/2018





