
## `ld-windowed`

This directory contains scripts used to calculate LD ($r^2$) across the mating type
locus. LD is then summarized in windows using the ZnS statistic (Kelly 1997).

### Transposing the aligned FASTA - `transpose_aligned_fasta.py`

For ease of LD calculations, we begin by transposing the aligned FASTA
from `alignment-lastz` once again. This script also computes whether a
site is usable for LD calculations and indicates that in an added column.

### $r^2$ across the mt locus - `r2_calc.py`

This script uses the transposed FASTA to calculate pairwise LD
for all SNPs in a given windowsize. 

### ZnS across the mt locus - `zns_calc.py`

Finally, this script sums the $r^2$ values generated above to 
calculate ZnS in non-overlapping windows.
