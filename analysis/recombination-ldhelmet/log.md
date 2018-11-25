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












