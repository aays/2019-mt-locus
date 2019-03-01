
## `ld-windowed`

This directory contains scripts used to calculate LD ($r^2$) across the mating type
locus. LD is then summarized in windows using the ZnS statistic (Kelly 1997).

### Transposing the aligned FASTA - `transpose_aligned_fasta.py`

For ease of LD calculations, we begin by transposing the aligned 'master' 
FASTA from `alignment-lastz` once again. This script also computes whether a
site is usable for LD calculations and indicates that in an added column.

### $r^2$ across the mt locus - `r2_calc.py`

This script uses the transposed FASTA to calculate pairwise LD
for all SNPs in a given windowsize. 

### ZnS across the mt locus - `zns_calc.py`

Finally, this script sums the $r^2$ values generated above to 
calculate ZnS in non-overlapping windows.

### `mt_locus_ld.sh`

This shell script runs the above pipeline on all alignment chunks
from the LASTZ alignment. Contains hardcoded paths.

### `r2_calc_region.py`

This is a quick and dirty modification of `r2_calc` that calculates r2 within a
given region. Useful for estimating r2 over long stretches of sequence (i.e. a
full chromosome) by estimating r2 for chunks in parallel.  Currently not very
efficient though - instead of using tabix, simply iterates through file until
the desired region is found...

### `zns_genes.py`

This script takes in a directory containing aligned FASTAs with individual
genes across all strains and performs the full LD pipeline (transposition,
r2 calculation, zns calculation). Creates a directory called `temp-transposed`
to store transposed files + r2 values in. This also requires a csv
(generated in `cds-popgen`) containing other gene-specific stats - the output
of this script is that same csv with a 'zns' column tacked on.
