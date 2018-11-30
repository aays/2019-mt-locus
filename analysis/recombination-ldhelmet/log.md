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














