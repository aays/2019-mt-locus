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

so:

`transpose_fastas.py` -
1. open transposed file to write to
2. for each fasta in `alignments`:
    a. for each position in the fasta:
        a. if mt+ is not '-': write out position and base for strains
        b. elif mt+ is '-': skip over this site

although to make sure we get the strain names,
we have to modify `align_mt_fasta_maf.py` to output
those as fasta ids (instead of just the regions)

```
mkdir -p data/aligned-fastas/alignments_strains
time python3.5 analysis/alignment-lastz/align_mt_fasta_maf.py \
--plus data/aligned-fastas/plus_strains_ref.fasta \
--minus data/aligned-fastas/minus_strains_ref.fasta \
--alignment data/alignment-lastz/lastz-align-10k-gapped.maf \
--bed data/alignment-lastz/lastz-align-10k-gapped-filtered.bed \
--outdir data/aligned-fastas/alignments_strains

time python3.5 analysis/alignment-lastz/transpose_fastas.py \
--directory data/aligned-fastas/alignments_strains/ \
--outfile data/aligned-fastas/mt_aligned_transposed.txt
```

`remove_duplicates.R` 
1. read in file above as df
2. dplyr::distinct() to remove duplicates
3. assert that no position appears twice

```bash
Rscript analysis/alignment-lastz/remove_duplicates.R \
data/aligned-fastas/mt_aligned_transposed.txt \
data/aligned-fastas/mt_aligned_transposed_filtered.txt
```

and then:

`combine_fastas.py` - this input file will have no mt+ gaps or duplicates
1. open transposed file
2. create dict with strain names from header
3. instantiate counter = 298299
3. for line in fasta:
    a. check that position = counter
    b. grow sequence string (dict value) for all strains
    c. if not position = counter # gap happened
        a. gap_size = position - counter
        b. increment sequence string with gap_size number of Ns
        c. then proceed to next iteration

```bash
time python3.5 analysis/alignment-lastz/combine_fastas.py \
--file data/aligned-fastas/mt_aligned_transposed_filtered.txt \
--outfile data/aligned-fastas/mt_aligned_final.fasta
```
