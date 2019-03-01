
## `organelle-linkage`

This directory contains scripts used to calculate LD (D') between mating
type loci and their corresponding organelle genomes under a uniparental
inheritance paradigm.

### Subsetting VCFs - `vcf_subset.py`

This general purpose script is used for subsetting VCFs based on either
or both of a list of chromosomes and a certain fraction of sites to keep.

### Inter-chromosomal LD - `ld_calc.py`

This script will calculate all three of the D, D', and r2 measures of LD
between all pairwise combinations of sites at two user-provided regions
in a single VCF file.

### LD between random pairs - `all_pairs_calc.py`

This script will calculate the LD stats listed above between _all_ pairs
of SNPs in a given VCF. 
