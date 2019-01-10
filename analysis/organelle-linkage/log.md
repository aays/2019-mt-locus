log for organelle linkage analysis

this has already been completed - just need to redo and put in this folder

## 5/1/2019

original shell script - https://github.com/aays/fiftyshadesofgreen/blob/master/vcf/correctstrains.sh

new GATK symlink:

```bash
ln -sv /opt/gatk/GenomeAnalysis.jar bin/
```

creating the necessary VCFs:

```bash
mkdir -p data/references/gvcfs/
mkdir -p data/organelle-linkage/vcfs/
cd data/references/gvcfs/
ln -sv /scratch/research/data/chlamydomonas/quebec/gVCFs/* .
ln -sv /scratch/research/data/chlamydomonas/species_wide/gVCFs/* .
cd ../../../

touch minus.intervals
for line in mtDNA mtMinus:286502-287123 mtMinus:308874-311600;
    do echo $line >> minus.intervals;
done

touch plus.intervals
for line in chromosome_6:480316-484847 chromosome_6:559654-561075 cpDNA;
    do echo $line >> plus.intervals;
done

java -jar ./bin/GenomeAnalysisTK.jar -T GenotypeGVCFs \
-R /scratch/research/references/chlamydomonas/5.3_chlamy_w_organelles_mt_minus/chlamy.5.3.w_organelles_mtMinus.fasta \
-L minus.intervals \
--variant data/references/gvcfs/CC2935.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC2938.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3059.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3061.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3062.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3063.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3073.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3075.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3079.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3084.haplotypeCalled.g.vcf.gz \
-o data/organelle-linkage/vcfs/mtmtd1midminus.vcf

java -jar ./bin/GenomeAnalysisTK.jar -T GenotypeGVCFs \
-R /scratch/research/references/chlamydomonas/5.3_chlamy_w_organelles_mt_minus/chlamy.5.3.w_organelles_mtMinus.fasta \
-L plus.intervals \
--variant data/references/gvcfs/CC2936.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC2937.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3060.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3064.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3065.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3068.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3071.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3076.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3086.haplotypeCalled.g.vcf.gz \
-o data/organelle-linkage/vcfs/cpmtafusplus.vcf

bgzip data/organelle-linkage/vcfs/mtmtd1midminus.vcf
tabix data/organelle-linkage/vcfs/mtmtd1midminus.vcf.gz

bgzip data/organelle-linkage/vcfs/cpmtafusplus.vcf
tabix data/organelle-linkage/vcfs/cpmtafusplus.vcf.gz    

rm plus.intervals
rm minus.intervals
```

now to migrate and use older scripts from the ld paper repo:

```bash
time python3.5 analysis/organelle-linkage/singlevcfcalc.py \
-v data/organelle-linkage/mtmtd1midminus.vcf.gz \
-r mtMinus mtDNA \
-l d/dprime/r2 > data/organelle-linkage/mtmtd1mid.txt

time python3.5 analysis/organelle-linkage/singlevcfcalc.py \
-v data/organelle-linkage/cpmtafusplus.vcf.gz \
-r chromosome_6 cpDNA \
-l d/dprime/r2 > data/organelle-linkage/cpmtafus.txt
``` 

so I'm seeing lots of D' = 0 in the mt- comparison, all of which is explained
by calls at CC3062, which are often at a far lower depth (~10-15) than the remaining
samples (DP = 30-80). going to try this again with a DP = 20 filter as well
in `popgen.straingetter`

update - this script is bonked in that if after filtering for strains a SNP is no
longer usable (ie becomes invariant) - straingetter can't handle it.

although frankly neither can I. time to start fresh with a new script:

```bash
time python3.5 analysis/organelle-linkage/ld_calc.py \
--vcf_file data/organelle-linkage/vcfs/mtmtd1midminus.vcf.gz \
--regions mtMinus mtDNA \
--outfile data/organelle-linkage/minus.txt

time python3.5 analysis/organelle-linkage/ld_calc.py \
--vcf_file data/organelle-linkage/vcfs/cpmtafusplus.vcf.gz \
--regions chromosome_6 cpDNA \
--outfile data/organelle-linkage/plus.txt
```

looks good!

mtDNA with itself:

```bash
time python3.5 analysis/organelle-linkage/ld_calc.py \
--vcf_file data/organelle-linkage/vcfs/mtmtd1midminus.vcf.gz \
--regions mtDNA mtDNA \
--outfile data/organelle-linkage/mt_only.txt
```

## 8/1/2019

now for the random inter-chr pairs calculation.

first, we'll make a vcf containing all the strains being used here:

```bash
touch exclusions.intervals
for line in mtDNA cpDNA mtMinus ;
    do echo $line >> exclusions.intervals;
done

java -jar ./bin/GenomeAnalysisTK.jar -T GenotypeGVCFs \
-R /scratch/research/references/chlamydomonas/5.3_chlamy_w_organelles_mt_minus/chlamy.5.3.w_organelles_mtMinus.fasta \
-XL exclusions.intervals \
--variant data/references/gvcfs/CC2935.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC2938.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3059.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3061.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3062.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3063.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3073.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3075.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3079.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3084.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC2936.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC2937.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3060.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3064.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3065.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3068.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3071.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3076.haplotypeCalled.g.vcf.gz \
--variant data/references/gvcfs/CC3086.haplotypeCalled.g.vcf.gz \
-o data/organelle-linkage/vcfs/all_variants.vcf
```

where `exclusions.intervals` is
```bash
cpDNA
mtDNA
mtMinus
```

finally done! (took 2.38 hours)

```bash
bgzip data/organelle-linkage/vcfs/all_variants.vcf
tabix data/organelle-linkage/vcfs/all_variants.vcf.gz
rm exclusions.intervals
```

now to subset:

```bash
time python3.5 analysis/organelle-linkage/vcf_subset.py \
--vcf data/organelle-linkage/vcfs/all_variants.vcf.gz \
--filter_fraction 0.0002 \
--outfile data/organelle-linkage/vcfs/all_variants_filtered.vcf

bgzip data/organelle-linkage/vcfs/all_variants_filtered.vcf
tabix data/organelle-linkage/vcfs/all_variants_filtered.vcf.gz
```

and, finally, the all pairs calculation:

```bash
time python3.5 analysis/organelle-linkage/all_pairs_calc.py \
--vcf_file data/organelle_linkage/vcfs/all_variants_filtered.vcf.gz \
--outfile data/organelle_linkage/allpairs.txt
```

aaand done!

