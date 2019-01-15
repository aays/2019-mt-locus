
## `recombination-ldhelmet`

This directory contains scripts used to calculate recombination rates
across the mt locus by applying LDhelmet (Chan et al, 2012) on FASTA
files generated in `alignment-lastz`.

### Running LDhelmet - `ldhelmet_indiv.sh`

This bash script is used to run the various steps of LDhelmet on
an individual FASTA file, such as the mt-aligned file previously generated.

### Correcting coordinates

Since the mt+ locus (which we are using as a 'reference track') starts
at position 298298 of chromosome 6, the coordinates in the LDhelmet outfile
need to be updated accordingly.

`ldhelmet_aligned_clean.py` will clean LDhelmet output for the aligned FASTA.

`ldhelmet_mt_only_clean.py` will clean LDhelmet output for the mt-separated FASTAs.

### Generating long-format recombination estimates - `generate_mt_long.py`

This script will combine the recombination rates calculated for the aligned mt loci
as well as those in the mt+ alone into a long-format file spanning the entire mt locus,
with individual positions as keys and per bp recombination rates as values.
Each position is also demarcated as being gametologous or non-gametologous.

This can be used for easy parsing and recombination rate lookup at any site
or region within the mt locus. 
