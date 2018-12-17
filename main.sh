# replicate analyses

# lastz alignment
echo "Aligning mt loci with LASTZ..."
time ./bin/lastz data/references/mt_plus.fasta data/references/mt_minus.fasta \
--output=data/alignment-lastz/lastz-align-10k-gapped.maf \
--hspthresh=10000 \
--format=maf

time ./bin/lastz data/references/mt_plus.fasta data/references/mt_minus.fasta \
--output=data/alignment-lastz/lastz-align-10k-gapped.bed \
--hspthresh=10000 \
--format=general
echo "Done."

# clean the bed file
Rscript analysis/alignment-lastz/clean_lastz_output.R \
-f data/alignment-lastz/lastz-align-10k-gapped.bed \
-o data/alignment-lastz/lastz-align-10k-gapped-filtered.bed

# aligning the fastas
echo "Creating mt-aligned fasta..."
mkdir -p data/aligned-fastas/alignments
time python3.5 analysis/alignment-lastz/align_mt_fasta_maf.py \
--plus data/aligned-fastas/plus_strains_ref.fasta \
--minus data/aligned-fastas/minus_strains_ref.fasta \
--alignment data/alignment-lastz/lastz-align-10k-gapped.maf \
--bed data/alignment-lastz/lastz-align-10k-gapped-filtered.bed \
--outdir data/aligned-fastas/alignments

echo "Creating mt-separated fastas..."
time python3.5 analysis/alignment-lastz/make_mt_only.py \
--fasta data/aligned-fastas/plus_strains_ref.fasta \
--alignment data/alignment-lastz/lastz-align-10k-gapped-filtered.bed \
--mt_allele plus \
--output data/aligned-fastas/plus_non_gametolog.fasta

time python3.5 analysis/alignment-lastz/make_mt_only.py \
--fasta data/aligned-fastas/minus_strains_ref.fasta \
--alignment data/alignment-lastz/lastz-align-10k-gapped-filtered.bed \
--mt_allele minus \
--output data/aligned-fastas/minus_non_gametolog.fasta


# recombination rate estimation
## gametolog recombination
echo "LDhelmet runs for aligned fasta..."
time bash analysis/recombination-ldhelmet/ldhelmet_init.sh
