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

## 8/1/2019

so a pretty serious error in the alignment has now been fixed (see
the `alignment-lastz` log - there should be a lot fewer mismatches
in the alignment now

redoing the above steps:

`r2_calc` took just 25 minutes, and `zns_calc` 9! much faster than the
multi-hour processes of the past - I guess that's what we get with
a bunch of false misalignments

now for LD in the plus non-gametolog file, to compare against:

```bash
time python3.5 analysis/ld-windowed/transpose_aligned_fasta.py \
--fasta data/aligned-fastas/plus_non_gametolog.fasta \
--outfile data/ld-windowed/plus_non_gametolog_long.txt \
--offset 298298

time python3.5 analysis/ld-windowed/r2_calc.py \
--filename data/ld-windowed/plus_non_gametolog_long.txt \
--windowsize 1000 \
--outfile data/ld-windowed/plus_NG_r2_1k.txt \
--non_gametolog

time python3.5 analysis/ld-windowed/zns_calc.py \
--filename data/ld-windowed/plus_NG_r2_1k.txt \
--windowsize 1000 \
--outfile data/ld-windowed/plus_NG_zns_1k.txt
```

## 11/1/2019

running the scripts above on the gametolog chunks individually:

```bash
mkdir -p data/ld-windowed/alignments-temp
mkdir -p data/ld-windowed/r2
mkdir -p data/ld-windowed/zns

for fname in data/aligned-fastas/alignments/*fasta; do
    base=$(basename $fname .fasta)

    time python3.5 analysis/ld-windowed/transpose_aligned_fasta.py \
    --fasta ${fname} \
    --outfile data/ld-windowed/alignments-temp/${base}_long.txt
    --offset 298298

    time python3.5 analysis/ld-windowed/r2_calc.py \
    --filename data/ld-windowed/alignments-temp/${base}_long.txt \
    --windowsize 1000 \
    --outfile data/ld-windowed/r2/${base}_r2_1k.txt

    time python3.5 analysis/ld-windowed/zns_calc.py \
    --filename data/ld-windowed/r2/${base}_r2_1k.txt \
    --windowsize 1000 \
    --outfile data/ld-windowed/zns/${base}_zns_1k.txt;
done

```

## 12/1/2019

so our transposition script won't work with these files.

two issues:
1. the offset needs to be dynamically updated from each file
2. these individual alignments still have plus gaps in them

perhaps the alignment script can be updated with an 
'infer-offset' setting? this setting would use the
(standardized) filenames to infer offset levels
and also crunch together gaps in the plus.

## 13/1/2019

trying this after updating the script:

```bash
time python3.5 analysis/ld-windowed/transpose_aligned_fasta.py \
--fasta data/aligned-fastas/alignments/chromosome_6_309976-323648.fasta \
--outfile test2.out \
--infer_offset
```

looks good!

round 2:

```bash
mkdir -p data/ld-windowed/alignments-temp
mkdir -p data/ld-windowed/r2
mkdir -p data/ld-windowed/zns

for fname in data/aligned-fastas/alignments/*fasta; do
    base=$(basename $fname .fasta)

    time python3.5 analysis/ld-windowed/transpose_aligned_fasta.py \
    --fasta ${fname} \
    --outfile data/ld-windowed/alignments-temp/${base}_long.txt
    --infer_offset

    time python3.5 analysis/ld-windowed/r2_calc.py \
    --filename data/ld-windowed/alignments-temp/${base}_long.txt \
    --windowsize 1000 \
    --outfile data/ld-windowed/r2/${base}_r2_1k.txt

done

```

we'll need to combine these r2 outfiles into a single zns calculation somehow:

```bash
time python3.5 analysis/ld-windowed/zns_calc.py \
--filename data/ld-windowed/r2/${base}_r2_1k.txt \
--windowsize 1000 \
--outfile data/ld-windowed/zns/${base}_zns_1k.txt;
```

update: there's an issue here, where overlapping alignments are returning different
LD values. if there are two alignments in the mt- plus orientation, they should have
the same level of LD, but a third alignment in the mt- minus orientation might result
in a different LD value.

we should resolve these as before - reducing multiple overlapping alignments
into the one alignment with the highest score.

```
> r2_files %>% filter(snp1 == 422115, snp2 == 422301)
# A tibble: 3 x 4
                                  name   snp1   snp2     r2
                                 <chr>  <int>  <int>  <dbl>
1 chromosome_6_420433-424443_r2_1k.txt 422115 422301 1.0000
2 chromosome_6_421359-424446_r2_1k.txt 422115 422301 0.7875
3 chromosome_6_421599-424445_r2_1k.txt 422115 422301 1.0000
```

back to `alignment-lastz`!


## 15/1/2019

while the LDhelmet script reruns, we could get started
on running the LD analysis across chromosome 6:

```bash
time python3.5 analysis/ld-windowed/transpose_aligned_fasta.py \
--fasta data/aligned-fastas/chromosome_6_all.fasta \
--outfile data/ld-windowed/chromosome_6_long.txt \
--offset 0

time python3.5 analysis/ld-windowed/r2_calc.py \
--filename data/ld-windowed/chromosome_6_long.txt \
--windowsize 1000 \
--outfile data/ld-windowed/r2/chromosome_6_r2_1k.txt

time python3.5 analysis/ld-windowed/zns_calc.py \
--filename data/ld-windowed/r2/chromosome_6_r2_1k.txt \
--windowsize 1000 \
--outfile data/ld-windowed/zns/chromosome_6_zns_1k.txt
```

## 16/1/2019

alright, this has been running for 24 hours now, and is
well past the first ~1m SNPs - going to stop this now,
since we only really need LD around the mt locus

now for zns:

```bash
time python3.5 analysis/ld-windowed/zns_calc.py \
--filename data/ld-windowed/r2/chromosome_6_r2_1k.txt \
--windowsize 1000 \
--outfile data/ld-windowed/zns/chromosome_6_zns_1k.txt
```

in the future, it'd be helpful to have a SAM-style
defined range in the r2 script (ie 1-100000)

we also need to combine all the individual zns files
into an 'mt locus-wide' file to plot with

```R
library(tidyverse)
library(fs)
library(magrittr)

fnames <- dir_ls('data/ld-windowed/zns', regexp = '.*[0-9]{6}.*')
zns_all <- map_dfr(fnames, read_delim, delim = ' ', col_types = cols())
# resolve cases of multiple values in the same window by
# selecting the window that has more sites
zns_all %<>%
    group_by(start) %>%
    filter(site_count == max(site_count)) %>%
    ungroup()
write_delim(zns_all, path = 'data/ld-windowed/mt_locus_zns.txt',
            delim = ' ')
```

or simply:

```bash
time Rscript analysis/ld-windowed/combine_zns.R \
--directory data/ld-windowed/zns \
--outfile data/ld-windowed/mt_locus_zns.txt
```

## 20/1/2019

need to make a zns script that reads in a data frame
of gene coordinates and returns zns for each

although it'd be possible to write a script that
looks up zns from existing r2 files, it's probably
a better idea to have this script do all
three steps at once for each set of coordinates
(i.e. transposition, r2, zns)

could help to import functions from the existing scripts?
we'll see I suppose

this script will need:
- input dataframe from `cds-popgen`
- dir containing reference fastas? or overall ref fasta (`translationaligner/curated`)
- position doesn't actually matter here - so no need for offset

no windowsize will be used here - ZnS is for all SNPs
within a given region

## 21/1/2019

let's try running the script:

```bash
time python3.5 analysis/ld-windowed/zns_genes.py \
--filename analysis/cds-popgen/polymorphism/polymorphism_all.csv \
--directory analysis/cds-popgen/TranslationAligner/curated/ \
--outfile test_zns.out
```

looks good! although it seems some genes didn't have
any variants at all... would be interesting to look at these
down the line.




