# replicate analyses

# add a log functionality later if possible
# ie echo 'Starting x' >> log.txt
# { time sleep 1 ; } 2>> log.txt
# this will keep track of times too
# touch log.txt

### LASTZ ALIGNMENT ###

echo "Aligning mt loci with LASTZ..."
time ./bin/lastz data/references/mt_plus.fasta data/references/mt_minus.fasta \
--output=data/alignment-lastz/lastz-align-30k-gapped.maf \
--hspthresh=30000 \
--format=maf

time ./bin/lastz data/references/mt_plus.fasta data/references/mt_minus.fasta \
--output=data/alignment-lastz/lastz-align-30k-gapped.bed \
--hspthresh=30000 \
--format=general
echo "Done."

# clean the bed file
time python3.5 analysis/alignment-lastz/clean_lastz_output.py \
--filename data/alignment-lastz/lastz-align-30k-gapped.bed \
--threshold 0.75 \
--outfile data/alignment-lastz/lastz-align-30k-gapped-filtered.bed

sleep 3

### CREATE MT LOCUS FASTAS ###

echo "Creating alignment files..."
mkdir -p data/aligned-fastas/alignments
time python3.5 analysis/alignment-lastz/align_mt_fasta_maf.py \
--plus data/aligned-fastas/plus_strains_ref.fasta \
--minus data/aligned-fastas/minus_strains_ref.fasta \
--alignment data/alignment-lastz/lastz-align-30k-gapped.maf \
--bed data/alignment-lastz/lastz-align-30k-gapped-filtered.bed \
--outdir data/aligned-fastas/alignments/

echo "Done."
echo "Transposing alignments..."
time python3.5 analysis/alignment-lastz/transpose_fastas.py \
--directory data/aligned-fastas/alignments/ \
--outfile data/aligned-fastas/mt_aligned_transposed.txt

Rscript analysis/alignment-lastz/remove_duplicates.R \
data/aligned-fastas/mt_aligned_transposed.txt \
data/aligned-fastas/mt_aligned_transposed_filtered.txt

echo "Done."
echo "Creating final mt-aligned file..."
time python3.5 analysis/alignment-lastz/combine_fastas.py \
--file data/aligned-fastas/mt_aligned_transposed_filtered.txt \
--outfile data/aligned-fastas/mt_aligned_final.fasta

echo "Creating C-domain appended file..."
time ./bin/vcf2fasta.py -v data/references/all_quebec.HC.vcf.gz \
-r data/references/mtPlus_ref.chromosome_6.fasta \
-i chromosome_6:826738-943474 \
-s CC2936 CC2937 CC3060 CC3064 CC3065 CC3068 CC3071 CC3076 CC3086 \
CC2935 CC2938 CC3059 CC3061 CC3062 CC3063 CC3073 CC3075 CC3079 CC3084 \
--min_GQ 30 > data/aligned-fastas/c_domain_aligned.fasta

cp data/aligned-fastas/c_domain_aligned.fasta \
data/aligned-fastas/alignments/chromosome_6_826738-943474.fasta

echo "Appending C domain..."
time python3.5 analysis/alignment-lastz/append_c_domain.py \
--mt_aligned data/aligned-fastas/mt_aligned_final.fasta \
--c_domain data/aligned-fastas/c_domain_aligned.fasta \
--outfile data/aligned-fastas/mt_aligned_all.fasta

echo "mt-aligned file created."
echo "Clearing intermediate files..."
rm -v data/aligned-fastas/mt_aligned_transposed*

sleep 3

# don't use the files below for rho estimation!
echo "Creating mt-separated fastas..."
time python3.5 analysis/alignment-lastz/make_mt_only.py \
--fasta data/aligned-fastas/plus_strains_ref.fasta \
--alignment data/alignment-lastz/lastz-align-30k-gapped-filtered.bed \
--mt_allele plus \
--output data/aligned-fastas/plus_non_gametolog.fasta

time python3.5 analysis/alignment-lastz/make_mt_only.py \
--fasta data/aligned-fastas/minus_strains_ref.fasta \
--alignment data/alignment-lastz/lastz-align-30k-gapped-filtered.bed \
--mt_allele minus \
--output data/aligned-fastas/minus_non_gametolog.fasta

# mask paralogous regions in plus ref
touch exclusions.intervals
for line in chromosome_6:535350-542150 chromosome_6:424800-427850; do
    echo ${line} >> exclusions.intervals;
done

time python3.5 analysis/alignment-lastz/mask_paralogs.py \
--filename data/aligned-fastas/plus_strains_ref.fasta \
--mask_intervals exclusions.intervals \
--offset 298298 \
--outfile data/aligned-fastas/plus_strains_ref_masked.fasta

# use masked file from here on out
mv -v data/aligned-fastas/plus_strains_ref.fasta data/aligned-fastas/plus_strains_ref_unmasked.fasta
mv -v data/aligned-fastas/plus_strains_ref_masked.fasta data/aligned-fastas/plus_strains_ref.fasta

sleep 3

### RHO ESTIMATION ###

echo "Done."
echo "Starting recombination rate estimation."

# recombination rate estimation
## gametolog recombination
echo "LDhelmet runs for aligned fasta..."
time bash analysis/recombination-ldhelmet/ldhelmet_indiv.sh \
data/aligned-fastas/mt_aligned_all.fasta

sleep 3

mv -v data/recombination-ldhelmet/recombination-estimates/mt_aligned_all.txt \
data/recombination-ldhelmet/recombination-estimates/mt_aligned_raw.txt # rename

echo "Correcting coordinates for LDhelmet output..."
time python3.5 analysis/recombination-ldhelmet/ldhelmet_aligned_clean.py \
-f data/recombination-ldhelmet/recombination-estimates/mt_aligned_raw.txt \
-o data/recombination-ldhelmet/recombination-estimates/mt_aligned_all.txt

## non-gametolog recombination
echo "LDhelmet runs for individual mt loci..."
time bash analysis/recombination-ldhelmet/ldhelmet_indiv.sh \
data/aligned-fastas/plus_strains_ref.fasta

sleep 3

time bash analysis/recombination-ldhelmet/ldhelmet_indiv.sh \
data/aligned-fastas/minus_strains_ref.fasta

sleep 3

for allele in plus minus; do
    time python3.5 analysis/recombination-ldhelmet/ldhelmet_mt_only_clean.py \
    --filename data/recombination-ldhelmet/recombination-estimates/${allele}_strains_ref.txt \
    --bed data/alignment-lastz/lastz-align-30k-gapped-filtered.bed \
    --allele ${allele} \
    --fasta data/aligned-fastas/${allele}_strains_ref.fasta \
    --outfile data/recombination-ldhelmet/recombination-estimates/${allele}_strains_ref_corrected.txt
done

echo "Done."

sleep 3

echo "Generating long-form recombination estimates..."
time python3.5 analysis/recombination-ldhelmet/generate_mt_long.py \
--mt_locus data/recombination-ldhelmet/recombination-estimates/mt_aligned_all.txt \
--plus data/recombination-ldhelmet/recombination-estimates/plus_strains_ref_corrected.txt \
--alignment data/alignment-lastz/lastz-align-30k-gapped-filtered.bed \
--fasta data/aligned-fastas/plus_strains_ref.fasta \
--outfile data/recombination-ldhelmet/recombination-estimates/mt_full_long.txt

echo "Done."

sleep 3

### MT LOCUS-WIDE LD ESTIMATIONS ###

echo "Begin LD calculation across mt locus."

mkdir -p data/ld-windowed/alignments-temp
mkdir -p data/ld-windowed/r2
mkdir -p data/ld-windowed/zns

count=1
filecount=$(ls data/aligned-fastas/alignments/*fasta | wc -l)

# loop over alignments
for fname in data/aligned-fastas/alignments/*fasta; do
    base=$(basename $fname .fasta)
    echo "Currently on ${base}"
    echo "File ${count} of ${filecount}"

    python3.5 analysis/ld-windowed/transpose_aligned_fasta.py \
    --fasta ${fname} \
    --outfile data/ld-windowed/alignments-temp/${base}_long.txt \
    --infer_offset

    sleep 1

    echo "Calculating r2..."
    python3.5 analysis/ld-windowed/r2_calc.py \
    --filename data/ld-windowed/alignments-temp/${base}_long.txt \
    --windowsize 1000 \
    --outfile data/ld-windowed/r2/${base}_r2_1k.txt

    sleep 1

    echo "Calculating ZnS..."
    python3.5 analysis/ld-windowed/zns_calc.py \
    --filename data/ld-windowed/r2/${base}_r2_1k.txt \
    --windowsize 1000 \
    --outfile data/ld-windowed/zns/chromosome_6_zns_1k.txt

    (( count ++ ))

done

rm data/ld-windowed/alignments-temp/*

sleep 3

### ORGANELLE LINKAGE ###

# mtDNA with mt-
time python3.5 analysis/organelle-linkage/ld_calc.py \
--vcf_file data/organelle-linkage/vcfs/mtmtd1midminus.vcf.gz \
--regions mtMinus mtDNA \
--outfile data/organelle-linkage/minus.txt

# cpDNA with mt+
time python3.5 analysis/organelle-linkage/ld_calc.py \
--vcf_file data/organelle-linkage/vcfs/cpmtafusplus.vcf.gz \
--regions chromosome_6 cpDNA \
--outfile data/organelle-linkage/plus.txt

# mtDNA with self
time python3.5 analysis/organelle-linkage/ld_calc.py \
--vcf_file data/organelle-linkage/vcfs/mtmtd1midminus.vcf.gz \
--regions mtDNA mtDNA \
--outfile data/organelle-linkage/mt_only.txt

# random inter-chr pairs
time python3.5 analysis/organelle-linkage/all_pairs_calc.py \
--vcf_file data/organelle_linkage/vcfs/all_variants_filtered.vcf.gz \
--outfile data/organelle_linkage/allpairs.txt

