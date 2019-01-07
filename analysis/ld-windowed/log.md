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

## 1/1/2019

took 11 hours - and it seems everything is in complete LD - this
script is bonked and needs work

fixed some bugs - attempt 2:

```bash
time python3.5 analysis/ld-windowed/r2_calc.py \
--filename data/ld-windowed/mt_aligned_long.txt \
--windowsize 25000 \
--outfile data/ld-windowed/mt_aligned_r2_25k.txt
```

another issue:
if there are mismatches in the alignment, such that
the long form file looks like this:

```
position S1 S2 S3 S4
50 A A T T
51 C C T T
52 G G A A
```

and so on, these will always be in full LD with
one another, because they are simply mismatches between
the alignment. the r2 script needs to be
rewritten to ignore these.

```bash
time python3.5 analysis/ld-windowed/r2_calc.py \
--filename data/ld-windowed/mt_aligned_long.txt \
--windowsize 20000 \
--outfile data/ld-windowed/mt_aligned_r2_20k.txt
```

## 2/1/2019

took nearly 10 hours but seems to have worked like a charm!

next up - need to write a script that calculates ZnS in non-overlapping windows


## 3/1/2019

the more I think about it, those pairs we were ignoring
should still be included. given

```
A A A T T T
C C C G G G
```

a signature of recombination would be something like

```
A T A T A T
C C C G G G
```

and so the first case still counts as usable data,
seeing as it indicates a lack of recombination.

`usable_pair` in the script should therefore
auto-return LD = 1 instead of claiming that the pair
is non-usable.

```bash
time python3.5 analysis/ld-windowed/r2_calc.py \
--filename data/ld-windowed/mt_aligned_long.txt \
--windowsize 1000 \
--outfile data/ld-windowed/mt_aligned_r2_1k.txt
```

took 4.65 hours

## 6/1/2019

alright - now for a script that calculates ZnS for a given windowsize,
given the output of `r2_calc.py`.

```bash
# test on smaller file
head -n 2001 data/ld-windowed/mt_aligned_r2_1k.txt > data/ld-windowed/mt_aligned_long_2000.txt

time python3.5 analysis/ld-windowed/zns_calc.py \
--filename data/ld-windowed/mt_aligned_long_2000.txt \
--windowsize 200 \
--outfile test_ld.txt
```

main command:

```bash
time python3.5 analysis/ld-windowed/zns_calc.py \
--filename data/ld-windowed/mt_aligned_r2_1k.txt \
--windowsize 1000 \
--outfile data/ld-windowed/mt_aligned_zns_1k.txt
```

wait - there's a serious bug in the r2 file:

```python
for record in [record1, record2]:
    plus_variants = list(set([record[key] for key in plus_strains]))
    minus_variants = list(set([record[key] for key in minus_strains]))
    if len(plus_variants) == 1 and len(minus_variants) == 1:
        return '1' # r2 = 1.0 by default
```

this means that even if just the first record's variation is
owing to a mismatch, LD between _both_ SNPs is automatically reported as 1.

running this again post-fix:

```bash
time python3.5 analysis/ld-windowed/r2_calc.py \
--filename data/ld-windowed/mt_aligned_long.txt \
--windowsize 1000 \
--outfile data/ld-windowed/mt_aligned_r2_1k_fix.txt
```

now took 5.5 hours - and now for the zns script:

```bash
time python3.5 analysis/ld-windowed/zns_calc.py \
--filename data/ld-windowed/mt_aligned_r2_1k.txt \
--windowsize 1000 \
--outfile data/ld-windowed/mt_aligned_zns_1k_fixed.txt
```

took 7.5 hours! 
