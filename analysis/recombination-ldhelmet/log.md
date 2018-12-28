log for LDhelmet recombination rate estimation

## 20/11/2018

now that the aligned fasta has been made in `analysis/alignment-lastz`, we can get to
using it as input for LDhelmet

the LDhelmet 1.9 executable already exists in `bin/` - the steps are:

1. `find_confs` - creates haplotype conf files
2. `table_gen` - creates lk lookup tables
3. `pade` - creates coefficient tables
4. `rjmcmc` - creates post files
5. `post_to_text` - creates tsv output

```bash
mkdir data/recombination-ldhelmet/intermediate-files
mkdir data/recombination-ldhelmet/recombination-estimates

time ./bin/ldhelmet find_confs \
--num_threads 10 \
--window_size 50 \
--output_file data/recombination-ldhelmet/intermediate-files/output.conf \
data/aligned-fastas/mt_aligned.fasta

time ./bin/ldhelmet table_gen \
--num_threads 10 \
--conf_file data/recombination-ldhelmet/intermediate-files/output.conf \
--theta 0.01 \
--rhos 0.0 0.1 10.0 1.0 100.0 \
--output_file data/recombination-ldhelmet/intermediate-files/output.lk

time ./bin/ldhelmet pade \
--num_threads 10 \
--conf_file data/recombination-ldhelmet/intermediate-files/output.conf \
--theta 0.01 \
--output_file data/recombination-ldhelmet/intermediate-files/output.pade

time ./bin/ldhelmet rjmcmc \
--num_threads 10 \
--window_size 50 \
--seq_file data/aligned-fastas/mt_aligned.fasta \
--lk_file data/recombination-ldhelmet/intermediate-files/output.lk \
--pade_file data/recombination-ldhelmet/intermediate-files/output.pade \
--num_iter 1000000 \
--burn_in 100000 \
--block_penalty 50 \
--output_file data/recombination-ldhelmet/intermediate-files/output.post

time ./bin/ldhelmet post_to_text \
--mean \
--perc 0.025 \
--perc 0.50 \
--perc 0.975 \
--output_file data/recombination-ldhelmet/recombination-estimates/mt_locus_recombination.txt
data/recombination-ldhelmet/intermediate-files/output.post

```
the above code has been saved to `ldhelmet_main.sh`

next up: LDhelmet on the mt+ and mt- individuals independently - `ldhelmet_alleles.sh`

## 21/11/2018

finally, we'll run LDhelmet on chromosome 6 independently:

```bash
time ./bin/vcf2fasta.py -v data/references/all_quebec.HC.vcf.gz \
-r data/references/mtMinus_ref.chromosome_6_and_mtMinus.fasta \
-i chromosome_6:1-9023763 \
-s CC2936 CC2937 CC3060 CC3064 CC3065 CC3068 CC3071 CC3076 CC3086 \
CC2935 CC2938 CC3059 CC3061 CC3062 CC3063 CC3073 CC3075 CC3079 CC3084 \
--min_GQ 30 > data/aligned-fastas/chromosome_6_all.fasta
```
followed by `ldhelmet_chr6.sh`

(I know I should ideally have a general LDhelmet script that takes in filenames at the command
line, instead of hardcoding paths in newer scripts... maybe something to retroactively make afterwards)

## 22/11/2018

added `mean_rho.R` - R script to parse all LDhelmet outputs in a directory
and return mean per bp rho estimates.

seems the aligned mt has the highest, while the mt alleles on their own
are comparable to chromosome 6's overall per bp recombination rate!

## 24/11/2018

redoing ldhelmet on the aligned fasta after correcting it -

(see log in analysis/alignment-lastz/)

first, I've removed all the previous ldhelmet intermediate files,
and renamed the final file to `mt_locus_old.txt`

and so:

```bash
time bash analysis/recombination-ldhelmet/ldhelmet_main.sh
```

# 25/11/2018

and here we are, back again:

```bash
time bash analysis/recombination-ldhelmet/ldhelmet_main.sh
time bash analysis/recombination-ldhelmet/ldhelmet_alleles.sh
```

## 29/11/2018

how do gapped characters in alignments affect LDhelmet?

```bash
cp mt_aligned.fasta mt_aligned_test.fasta
```

in the new copy, I've manually added several gap characters in the first sequence:

```
1 >CC3071|chromosome_6:298299-826737
2 ---------------------------------------------------------------------GGCGCTCACGGGGGTGCTGCCCTTCAGGGGCGCCAA
3 >CC3076|chromosome_6:298299-826737
4 GATTACTGTGGGCCGCGAGTACGATGGCACCAGTGTGGACATCTGGTCCATGGGCGTCATCCTGTACGAGGCGCTCACGGGGGTGCTGCCCTTCAGGGGCGCCAA
```

running ldhelmet using all the steps from `ldhelmet_main.sh` on this new file,
I get a substantial difference in recombination estimates:

```R
> map(d, weighted_mean)
$d_new
# A tibble: 1 x 3
  length weighted_rho  rho_per_bp
   <int>        <dbl>       <dbl>
1 528249     651.6867 0.001233673

$d_old
# A tibble: 1 x 3
  length weighted_rho   rho_per_bp
   <int>        <dbl>        <dbl>
1 528249     166.4918 0.0003151767
```

last test, just to be absolutely sure - how does having Ns instead of dashes affect things? 
(vcf2fasta generates Ns, not dashes)

```bash
cp mt_aligned_test.fasta mt_aligned_Ns.fasta
```

in this file, the exact same dashes are replaced with Ns.

after running ldhelmet on this again, seems we get the same result:

```R
> library(tidyverse)
Loading tidyverse: ggplot2
Loading tidyverse: tibble
Loading tidyverse: tidyr
Loading tidyverse: readr
Loading tidyverse: purrr
Loading tidyverse: dplyr
Conflicts with tidy packages ---------------------------------------------------
filter(): dplyr, stats
lag():    dplyr, stats
> ldhelmet_cols <- c('left_snp', 'right_snp', 'mean', 'p025', 'p50', 'p975')
> d_new_dash <- read_delim('ldhelmet-gap-temp/output.txt', delim = ' ', skip = 3, col_types = cols(), col_names = ldhe+ )t_cols()
Error in ldhelmet_cols() : could not find function "ldhelmet_cols"
> d_new_dash <- read_delim('ldhelmet-gap-temp/output.txt', delim = ' ', skip = 3,
+ col_types = cols(), col_names = ldhelmet_cols)
> d_new_N <- read_delim('ldhelmet-gap-temp/outputN.txt', delim = ' ', skip = 3,
+ col_types = cols(), col_names = ldhelmet_cols)
> d_old <- read_delim('data/recombination-ldhelmet/recombination-estimates/mt_locus_recombination.txt', delim = ' ',
+ skip = 3, col_types = cols(), col_names = ldhelmet_cols)
> weighted_mean <- function(df) {
+     output <- df %>%
+         transmute(length = right_snp - left_snp,
+                   weighted_rho = mean * length) %>%
+         summarise_all(sum) %>%
+         mutate(rho_per_bp = weighted_rho / length)
+     return(output)
+ }
> d <- list(d_new_dash, d_new_N, d_old)
> names(d) <- c('d_new_dash', 'd_new_N', 'd_old')
> map(d, weighted_mean)
$d_new_dash
# A tibble: 1 x 3
  length weighted_rho  rho_per_bp
   <int>        <dbl>       <dbl>
1 528249     651.6867 0.001233673

$d_new_N
# A tibble: 1 x 3
  length weighted_rho  rho_per_bp
   <int>        <dbl>       <dbl>
1 528249     651.6867 0.001233673

$d_old
# A tibble: 1 x 3
  length weighted_rho   rho_per_bp
   <int>        <dbl>        <dbl>
1 528249     166.4918 0.0003151767
```

so Ns should be fine when we're making the aligned fasta file. 


## 30/11/2018

so now that the mt separated files have been fixed, we can
get to running ldhelmet on both of them.

running `ldhelmet_alleles.sh` after editing:

```bash
time bash analysis/recombination-ldhelmet/ldhelmet_alleles.sh
```

## 16/12/2018

starting over - we need new scripts that iterate through the contents of
`data/aligned-fastas/alignments`, and do the following:

`ldhelmet_init.py` - for each file:
1. create a temp dir for ldhelmet run, and a final dir for outfiles
2. run five ldhelmet steps on input fasta
3. delete temp files, preserving only outfile

`ldhelmet_clean.py` - for each outfile:
1. parse through file
2. modify coordinates to match actual mt+ positions
3. rewrite output file with new coordinates to new outfile-only dir

`ldhelmet_combine.R` - for each outfile:
1. combine all files in new dir into single df
2. remove all duplicates
3. write final full-mt outfile

we also need scripts that do this for the mt-separated files:
1. run LDhelmet on entire masked mt locus
2. read in values and return mean rho value ignoring the masked regions

first attempt at new script:

```bash
time python3.5 analysis/recombination-ldhelmet/ldhelmet_init.py \
--directory data/aligned-fastas/alignments \
--ldhelmet bin/ldhelmet \
--outdir data/recombination-ldhelmet/recombination-estimates/
```

this seems to be clogging up memory - let's try 
`subprocess.call` (line 100)

if this works, push to repo

update: it didn't work - write a shell script tomorrow
instead that uses basename etc to do the same thing

make sure it also clears out temp files!

## 17/12/2018

first attempt at shell script equivalent:

```bash
time bash analysis/recombination-ldhelmet/ldhelmet_init.sh
```

this seems to hang up after 6 or so files.

let's individually run the one it messed up on.
we'll start by making a script that does the same
thing, but for only one input file

```bash
time bash analysis/recombination-ldhelmet/ldhelmet_indiv.sh \
data/aligned-fastas/alignments/chromosome_6_334439-334705.fasta
```

just froze on find_confs forever? what the shit?

in the meantime, the other script froze on a short 
alignment at 320000. although the alignment is ~200 bp,
the timestamp on the post file indicates LDhelmet just
stopped writing to the file at some point and sat there,
suddenly consuming 30 or so jobs on the server
instead of the 1-2 typical of an rjmcmc with this
many threads specified.

let's try this alignment individually and see if that makes
a difference:

```bash
time bash analysis/recombination-ldhelmet/ldhelmet-indiv.sh \
data/aligned-fastas/alignments/chromosome_6_320453-320614.fasta
```

if this works, we'll have to create a list of files in the
alignments directory and then use `split` to stick it into
chunks of ~3 filenames at a time to then run a modified LDhelmet
script on using the `while read -r` syntax in bash. although that'd
suck and mean more work, it also means the task can be parallelized
somewhat, which should save time in the long run.

update: it froze again. dammit. maybe LDhelmet can't
handle files that are this short? this one hardly spans
200 bp, and has only 24 SNPs. 

let's try to find an even shorter alignment and see if
the same issue repeats itself:

```bash
# 166 bp
time bash analysis/recombination-ldhelmet/ldhelmet-indiv.sh \
data/aligned-fastas/alignments/chromosome_6_374921-375087.fasta
```

although, looking back at the timestamps of the files from yesterday,
they got past this file...

```bash
$ ll
total 404K
-rw-r--r-- 1 hasans11 researchers  34K Dec 17 16:12 chromosome_6_298298-310006.txt
-rw-r--r-- 1 hasans11 researchers  53K Dec 17 16:26 chromosome_6_309976-323648.txt
-rw-r--r-- 1 hasans11 researchers 120K Dec 17 17:02 chromosome_6_310768-317056.txt
-rw-r--r-- 1 hasans11 researchers  38K Dec 17 17:19 chromosome_6_319517-320643.txt
-rw-r--r-- 1 hasans11 researchers 1.5K Dec 17 12:29 chromosome_6_320453-320614.txt
-rw-r--r-- 1 hasans11 researchers  28K Dec 17 12:41 chromosome_6_323640-330370.txt
-rw-r--r-- 1 hasans11 researchers  43K Dec 17 12:56 chromosome_6_325814-328093.txt
-rw-r--r-- 1 hasans11 researchers  70K Dec 17 13:13 chromosome_6_331467-348681.txt
```
notice how the 320k file took 17 minutes, but finished eventually.

let's leave this individual file that's currently underway
running on the server. if it eventually wraps up, 
we'll find a way to parallelize this as mentioned above,
since the only real problem is that it's being very slow. 

rjmcmc started at ~17:59

update: it's now 21:00 and it's still not done. something
is clearly very, very wrong.

on the flip side, let's try an LDhelmet run for a region
that's at least 10 kb - does this work while the shorter
files don't? could this have to do with how LDhelmet
calculates these recombination rates? for context - the 
example fasta bundled w/ LDhelmet is 25 kb

also changed all max threads back to 10 to prevent
the server from being overloaded

```bash
# started at 21:18
time bash analysis/recombination-ldhelmet/ldhelmet_indiv.sh \
data/aligned-fastas/alignments/chromosome_6_298298-310006.fasta
```

also, tomorrow, check whether the coordinates in the
outfiles are including or ignoring gaps! if the gaps
are being included, we'll have to update the `ldhelmet_clean.py`
script to test for that

also try running the smaller files with LDhelmet 1.7 instead?

## 18/12/2018

so this was done in ~11 minutes, which tells me
LDhelmet can't quite handle shorter sequences.

also, it appears the LDhelmet coordinates include
gaps - these need to be removed in the `ldhelmet_clean.py` script.

...you might have to remake the fastas, using a modified
`align_mt_fasta.maf` script that basically jumps over gaps
in the mt+.

an easier method would be to 'transpose' the fastas,
using the fastas that exist in `alignments` - these could
be written out like so:

```bash
position CCXXXX CCXXXX CCXXXX etc
298299 A A T
298300 G G G
```

and so on. this could later be read into R
for the removal of duplicates, after which we
could transpose them back out to a single fasta.

the method for this is detailed in the `alignment-lastz` log

anyways - IT WORKED! thank you sleepy ahmed for the epiphany

ldhelmet time:

```bash
time bash analysis/recombination-ldhelmet/ldhelmet_indiv.sh \
data/aligned-fastas/mt_aligned_final.fasta
```

and now to correct coordinates:

```bash
mv -v data/recombination-estimates/mt_aligned_final.txt \
data/recombination-estimates/mt_aligned_raw.txt

time python3.5 analysis/recombination-ldhelmet/ldhelmet_overall_clean.py \
-f data/recombination-ldhelmet/recombination-estimates/mt_aligned_raw.txt \
-o data/recombination-ldhelmet/recombination estimates/mt_aligned_final.txt
```

next, run LDhelmet on the two mt-separated files individually:

```bash
time bash analysis/recombination-ldhelmet/ldhelmet_indiv.sh \
data/aligned-fastas/plus_non_gametolog.fasta

time bash analysis/recombination-ldhelmet/ldhelmet_indiv.sh \
data/aligned-fastas/minus_non_gametolog.fasta
```

## 19/12/2018

next, we have to filter the mt-separated LDhelmet outfiles
for gametologous regions. sometimes, LDhelmet will bizarrely return
a rho value for a giant gametologous region of all Ns.

`ldhelmet_mt_only_clean.py`
1. takes in filtered bed file and LDhelmet outfile (+fasta)
2. use sets and range objects to reduce all bed intervals to a list
3. for each line in LDhelmet outfile
    if left_snp -> right_snp boundary covers value in bed intervals
        skip over the line
    else
        write line to new filtered outfile

```bash
time python3.5 analysis/recombination-ldhelmet/ldhelmet_mt_only_clean.py \
--filename data/recombination-ldhelmet/recombination-estimates/plus_non_gametolog.txt \
--bed data/alignment-lastz/lastz-align-10k-gapped-filtered.bed \
--allele plus \
--fasta data/aligned-fastas/plus_strains_ref.fasta \
--outfile data/recombination-ldhelmet/recombination-estimates/plus_non_gametolog_filtered.txt

time python3.5 analysis/recombination-ldhelmet/ldhelmet_mt_only_clean.py \
--filename data/recombination-ldhelmet/recombination-estimates/minus_non_gametolog.txt \
--bed data/alignment-lastz/lastz-align-10k-gapped-filtered.bed \
--allele minus \
--fasta data/aligned-fastas/minus_strains_ref.fasta \
--outfile data/recombination-ldhelmet/recombination-estimates/minus_non_gametolog_filtered.txt
```

and now to get mean rho:

```
Rscript analysis/recombination-ldhelmet/mean_rho.R \
-d data/recombination-ldhelmet/recombination-estimates/ \
-o data/recombination-ldhelmet/mean_rho.txt
```

wait - this increased the recombination estimates for the non-gametolog regions??

it seems a lot of this (at least in the minus) is being driven by a few regions:

```R
# A tibble: 5,236 x 8
   left_snp right_snp     mean   p0.025   p0.500   p0.975 length weighted_rho
      <int>     <int>    <dbl>    <dbl>    <dbl>    <dbl>  <int>        <dbl>
 1    51880     56974 2.158200 0.374070 1.962800 5.036900   5094  10993.87080
 2   204689    204843 2.947400 0.433570 2.753100 6.272700    154    453.89960
 3    51767     51880 2.158200 0.374070 1.962800 5.036900    113    243.87660
 4   105936    106081 1.667400 0.395820 1.502500 4.272500    145    241.77300
 5    57046     57154 2.158200 0.374070 1.962800 5.036900    108    233.08560
 6    51643     51731 2.158200 0.374070 1.962800 5.036900     88    189.92160
 7    57459     57535 2.158200 0.374070 1.962800 5.036900     76    164.02320
 8    57280     57341 2.158200 0.374070 1.962800 5.036900     61    131.65020
 9    57237     57280 2.158200 0.374070 1.962800 5.036900     43     92.80260
10    27520     31230 0.023593 0.013215 0.022049 0.042326   3710     87.53003
# ... with 5,226 more rows
```

what if we use a higher block penalty?

```bash
cd data/recombination-ldhelmet/recombination-estimates
mv minus_non_gametolog.txt minus_non_gametolog_init.txt

time bash analysis/recombination-ldhelmet/ldhelmet_indiv.sh \
data/aligned-fastas/minus_non_gametolog.fasta
```

in the meantime, looks like the python script generating non-gametolog
regions is flawed somehow, since this region is included in the alignment:

```R
> library(readr); library(magrittr); library(dplyr, warn.conflicts = FALSE)
> d <- read_tsv('../alignment-lastz/lastz-align-10k-gapped-filtered.bed', col_types = cols())
> d %<>% arrange(zstart2)
> d %<>% select(zstart2, end2, strand2)
> d %>% filter(zstart2 > 50000, zstart2 < 60000)
# A tibble: 2 x 3
  zstart2  end2 strand2
    <int> <int>   <chr>
1   52423 56995       +
2   58298 80105       +
```

so the higher block penalty did reduce rho:

```R
> d %>% summarise(rho = sum(weighted), length = sum(length)) %>%
+ mutate(out = rho / length)
# A tibble: 1 x 3
       rho length         out
     <dbl>  <int>       <dbl>
1 1001.304 321726 0.003112289
```

and annihilate that problematic region (which shouldn't be there!)

```R
> d %>% filter(left_snp == 51880)
# A tibble: 1 x 8
  left_snp right_snp       mean       p025        p50       p975 length
     <int>     <int>      <dbl>      <dbl>      <dbl>      <dbl>  <int>
1    51880     56974 9.2477e-06 6.3515e-07 7.7162e-06 1.9274e-05   5094
# ... with 1 more variables: weighted <dbl>
```

but we still need to correct the ldhelmet filtering script to
remove this region entirely.

testing the function manually gives the expected output:

```python
>>> current_interval = set(list(range(51880 - 1, 56974)))
>>> current_interval.issubset(intervals)
False
```

somehow just running the script again made the problematic values disappear?
but the recombination rate in that file is still stupidly high, while the denominator
for the aligned mt is wrong (spans the whole mt locus and not just gametolog length)

might just have to make one of those `is_gametolog` type files that can later be
read into an R data frame. 

```bash
time python3.5 analysis/recombination-ldhelmet/generate_mt_long.py \
--mt_locus data/recombination-ldhelmet/recombination-estimates/mt_aligned_final.txt \
--plus data/recombination-ldhelmet/recombination-estimates/plus_non_gametolog_filtered.txt \
--alignment data/alignment-lastz/lastz-align-10k-gapped-filtered.bed \
--fasta data/aligned-fastas/plus_strains_ref.fasta \
--outfile test_long.txt
```

so this script worked (and will be useful for plotting) but look:

```R
> d %>% group_by(is_gametolog) %>% summarise(mean_rho = mean(rho, na.rm = TRUE))
# A tibble: 2 x 2
  is_gametolog    mean_rho
         <int>       <dbl>
1            0 0.011467309
2            1 0.001986462
```

that doesn't make sense at all - something's off here

wait - instead of doing these on the masked files, why not run LDhelmet on the full
unmasked mt-only files and then filtering out parts that are covered by the lastz alignment?

```bash
time bash analysis/recombination-ldhelmet/ldhelmet_indiv.sh \
data/aligned-fastas/plus_strains_ref.fasta

time bash analysis/recombination-ldhelmet/ldhelmet_indiv.sh \
data/aligned-fastas/minus_strains_ref.fasta
```

first, we have to update the coordinates with a new script:

```bash
time python3.5 analysis/recombination-ldhelmet/ldhelmet_mt_full_clean.py \
--filename data/recombination-ldhelmet/recombination-estimates/plus_strains_ref.txt \
--allele plus \
--outfile data/recombination-ldhelmet/recombination-estimates/plus_strains_ref_corrected.txt
```

and so:

```bash
time python3.5 analysis/recombination-ldhelmet/generate_mt_long.py \
--mt_locus data/recombination-ldhelmet/recombination-estimates/mt_aligned_final.txt \
--plus data/recombination-ldhelmet/recombination-estimates/plus_strains_ref_corrected.txt \
--alignment data/alignment-lastz/lastz-align-10k-gapped-filtered.bed \
--fasta data/aligned-fastas/plus_strains_ref.fasta \
--outfile test_long_full.txt
```

this sort of remedied the problem - but only sort of 

```R
> d %>% group_by(is_gametolog) %>% summarise(m = mean(rho, na.rm = TRUE))
# A tibble: 2 x 2
  is_gametolog           m
         <int>       <dbl>
1            0 0.005526570
2            1 0.001986462
>
```

but this different value for `is_gametolog == 0` goes to show that
removing the masking actually does make a difference. 

it seems this huge jump is being driven by just two regions:

```R
> d %>% filter(is_gametolog == 0) %>% arrange(-rho) %>% group_by(rho) %>% tally() %>% arrange(-rho)
# A tibble: 96 x 2
        rho     n
      <dbl> <int>
 1 0.099478  6753
 2 0.077931  2135
 3 0.059614   270
 4 0.054090    95
 5 0.042391     5
 6 0.040939     7
 7 0.031879    18
 8 0.026011    20
 9 0.024268    18
10 0.020789   224
```

removing these makes a substantial difference -

```R
> d %>% filter(rho != 0.099478, rho != 0.077931) %>% group_by(is_gametolog) %>%
+ summarise(m = mean(rho, na.rm = TRUE))
# A tibble: 2 x 2
  is_gametolog           m
         <int>       <dbl>
1            0 0.000622426
2            1 0.001986462
```

perhaps we should redo *all* of this with a block penalty of 100 and see if that makes a difference.
although this will also necessitate backing it up somehow?

the sequence itself (for the first large chunk) looks largely uniform and has a low SNP density

```python
>>> show_aln(535380, 535450)
CC2936|chromosome_6:298299-826737 TCCAAGCACTGTTGCAGTTCCCTAGCGGGCCGGCCCAGCACTGCTGGATGCCTGTGACATGCGCCACAGC
CC2937|chromosome_6:298299-826737 TCCAAGCACTGTTGCAGTTCCCTAGCGGGCCGGCCCAGCACTGCTGGATGCCTGTGACATGCGCCACAGC
CC3060|chromosome_6:298299-826737 TCCAAGCACTGTTGCAGTTCCCTAGCGGGCCGGCCCAGCACTGCTGGATGCCTGTGACATGCGCCACAGC
CC3064|chromosome_6:298299-826737 TCCAAGCACTGTTGCAGTTCCCTAGCGGGCCGGCCCAGCACTGCTGGATGCCTGTGACATGCGCCACAGC
CC3065|chromosome_6:298299-826737 TCCAAGCACTGTTGCAGTTCCCTAGCGGGCCGGCCCAGCACTGCTGGATGCCTGTGACATGCGCCACAGC
CC3068|chromosome_6:298299-826737 TCCAAGCTCTGTTGCAGTTCCCTAGCGGGCCGGCCCAGCACTGCTGGATGCCTGTGACATGCGCCACAGC
CC3071|chromosome_6:298299-826737 TCCAAGCACTGTTGCAGTTCCCTAGCGGGCCGGCCCAGCACTGCTGGATGCCTGTGACATGCGCCACAGC
CC3076|chromosome_6:298299-826737 TCCAAGCACTGTTGCAGTTCCCTAGCGGGCCGGCCCAGCACTGCTGGATGCCTGTGACATGCGCCACAGC
CC3086|chromosome_6:298299-826737 TCCAAGCTCTGTTGCAGTTCCCTAGCGGGCCGGCCCAGCACTGCTGGATGCCTGTGACATGCGCCACAGC
```

we'll use a new script called `ldhelmet_indiv_100.sh`:

```bash
for filename in mt_aligned_final plus_strains_ref minus_strains_ref; do
    time bash analysis/recombination-ldhelmet/ldhelmet_indiv_100.sh \
    data/aligned-fastas/${filename}.fasta;
done
```

## 20/12/2018

today:
- check on LDhelmet results from yesterday
- run LASTZ with the mt- as the target - does this cover the high rho plus parts?

### LDhelmet results

let's correct the coordinates of our plus outfile:


```bash
mv data/recombination-ldhelmet/recombination-estimates/mt_aligned_final_100.txt \
data/recombination-ldhelmet/recombination-estimates/mt_aligned_raw_100.txt 

time python3.5 analysis/recombination-ldhelmet/ldhelmet_mt_full_clean.py \
--filename data/recombination-ldhelmet/recombination-estimates/plus_strains_ref_100.txt \
--allele plus \
--outfile data/recombination-ldhelmet/recombination-estimates/plus_strains_ref_100_corrected.txt

time python3.5 analysis/recombination-ldhelmet/ldhelmet_overall_clean.py \
-f data/recombination-ldhelmet/recombination-estimates/mt_aligned_raw.txt \
-o data/recombination-ldhelmet/recombination-estimates/mt_aligned_final.txt

time python3.5 analysis/recombination-ldhelmet/generate_mt_long.py \
--mt_locus data/recombination-ldhelmet/recombination-estimates/mt_aligned_final_100.txt \
--plus data/recombination-ldhelmet/recombination-estimates/plus_strains_ref_100_corrected.txt \
--alignment data/alignment-lastz/lastz-align-10k-gapped-filtered.bed \
--fasta data/aligned-fastas/plus_strains_ref.fasta \
--outfile test_long_full_100.txt
```

seems LDhelmet flattened out most of the variation in rho:

```R
> d %>% group_by(is_gametolog) %>%
+ summarise(weight = mean(rho, na.rm = TRUE))
# A tibble: 2 x 2
  is_gametolog      weight
         <int>       <dbl>
1            0 0.001335540
2            1 0.001522998
```

I don't think block = 100 is the way to go after all. 

in `alignment-lastz` we'll now try the mt- target alignment

## 24/12/2018

so we've had a bit of an overhaul (see `alignment-lastz` log)
where we've switched to a different lastz score threshold,
altering our original alignment dataset. earlier analyses
could (fortunately) be easily rerun with `main.sh` - once
this is done, we can use the long format files to precisely
calculate recombination rates in and out of gametologs. 

there are two versions of this going for unique regions - one 
where non-gametolog recombination is calculated from the full
fastas, while the other has recombination calculated from fastas
where the gametologs are masked. I expect rho to (likely falsely)
be higher in the second, since in a region with simply less
data, per base recombination rates will be inflated. 

## 25/12/2018

merry christmas!

it seems we have two regions inflating the crap out of
the non-gametolog regions.

```R
d_ref %>%
+ filter(is_gametolog == 0) %>%
+ group_by(rho) %>%
+ tally() %>%
+ arrange(desc(rho))
# A tibble: 103 x 2
        rho     n
      <dbl> <int>
 1 0.099478  6753
 2 0.077931  2979
 3 0.059614   270
 4 0.054090    95
 5 0.042391     5
 6 0.040939     7
 7 0.031879    18
 8 0.026011    20
 9 0.024268    18
10 0.020789   224
# ... with 93 more rows

> d_ref %>% filter(rho != 0.099478, rho != 0.077931) %>%
+ group_by(is_gametolog) %>%
+ summarise(mean_rho = mean(rho, na.rm = TRUE))
# A tibble: 2 x 2
  is_gametolog     mean_rho
         <int>        <dbl>
1            0 0.0005454173
2            1 0.0014719018
> d_ref %>%
+ group_by(is_gametolog) %>%
+ summarise(mean_rho = mean(rho, na.rm = TRUE))
# A tibble: 2 x 2
  is_gametolog    mean_rho
         <int>       <dbl>
1            0 0.003787625
2            1 0.001471902
```

these regions are 535387-542139 and 424836-427814 in
the mt+ locus - when eyeballing the GFF, there aren't
any genes in either of these regions.

what does base coverage look like here? back to
the `alignment-lastz` log

## 26/12/2018

we've re-redone the recombination calcalation, this time
actually using the correct VCF:

but the same two regions again:

```R
> d %>%
+ group_by(is_gametolog) %>%
+ summarise(mean_rho = mean(rho, na.rm = TRUE))
# A tibble: 2 x 2
  is_gametolog    mean_rho
         <int>       <dbl>
1            0 0.008301928
2            1 0.002693016

> d %>% filter(is_gametolog == 0) %>%
+ group_by(rho) %>% tally() %>% arrange(desc(rho))
# A tibble: 157 x 2
        rho     n
      <dbl> <int>
 1 0.292340  6753
 2 0.117260   169
 3 0.090423    37
 4 0.081629     3
 5 0.049919  3298
 6 0.029169   614
 7 0.027580    35
 8 0.027531     7
 9 0.027299   267
10 0.026506     4
# ... with 147 more rows

> d %>% filter(rho != 0.292340, rho != 0.049919) %>%
+ group_by(is_gametolog) %>%
+ summarise(mean_rho = mean(rho, na.rm = TRUE))
# A tibble: 2 x 2
  is_gametolog     mean_rho
         <int>        <dbl>
1            0 0.0006071983
2            1 0.0026930156
```

what the hell is going on in these regions? back
to the `alignment-lastz` log





