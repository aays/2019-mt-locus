# replicate analyses

# add a log functionality later if possible
# ie echo 'Starting x' >> log.txt
# { time sleep 1 ; } 2>> log.txt
# this will keep track of times too
# touch log.txt

# lastz alignment
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
Rscript analysis/alignment-lastz/clean_lastz_output.R \
-f data/alignment-lastz/lastz-align-30k-gapped.bed \
-o data/alignment-lastz/lastz-align-30k-gapped-filtered.bed

sleep 3

# aligning the fastas
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

sleep 3

echo "Done."
echo "Starting recombination rate estimation."

# recombination rate estimation
## gametolog recombination
echo "LDhelmet runs for aligned fasta..."
time bash analysis/recombination-ldhelmet/ldhelmet_indiv.sh \
data/aligned-fastas/mt_aligned_final.fasta

sleep 3

mv -v data/recombination-ldhelmet/recombination-estimates/mt_aligned_final.txt \
data/recombination-ldhelmet/recombination-estimates/mt_aligned_raw.txt # rename

echo "Correcting coordinates for LDhelmet output..."
time python3.5 analysis/recombination-ldhelmet/ldhelmet_aligned_clean.py \
-f data/recombination-ldhelmet/recombination-estimates/mt_aligned_raw.txt \
-o data/recombination-ldhelmet/recombination-estimates/mt_aligned_final.txt

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
--mt_locus data/recombination-ldhelmet/recombination-estimates/mt_aligned_final.txt \
--plus data/recombination-ldhelmet/recombination-estimates/plus_strains_ref_corrected.txt \
--alignment data/alignment-lastz/lastz-align-30k-gapped-filtered.bed \
--fasta data/aligned-fastas/plus_strains_ref.fasta \
--outfile data/recombination-ldhelmet/recombination-estimates/mt_full_long.txt

echo "Done."

sleep 3

echo "Begin LD calculation across mt locus."
echo "Transposing aligned fasta...'
time python3.5 analysis/ld-windowed/transpose_aligned_fasta.py \
--fasta data/aligned-fastas/mt_aligned_final.fasta \
--outfile data/ld-windowed/mt_aligned_long.txt \
--offset 298298 # mt plus coordinates

echo "Done."

sleep 3

echo "Calculating pairwise LD in 1 kb windows..."
time python3.5 analysis/ld-windowed/r2_calc.py \
--filename data/ld-windowed/mt_aligned_long.txt \
--windowsize 1000 \
--outfile data/ld-windowed/mt_aligned_r2_1k.txt

echo "Done."

sleep 3

echo "Computing ZnS statistic in 1 kb windows..."
time python3.5 analysis/ld-windowed/zns_calc.py \
--filename data/ld-windowed/mt_aligned_r2_1k.txt \
--windowsize 1000 \
--outfile data/ld-windowed/mt_aligned_zns_1k.txt
echo "ZnS calculation completed."
