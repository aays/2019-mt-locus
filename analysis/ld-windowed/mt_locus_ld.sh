mkdir -p data/ld-windowed/alignments-temp
mkdir -p data/ld-windowed/r2
mkdir -p data/ld-windowed/zns

count=1
filecount=$(ls data/aligned-fastas/alignments/*fasta | wc -l)

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
    --outfile data/ld-windowed/zns/${base}_zns_1k.txt

    (( count ++ ))

done

rm data/ld-windowed/alignments-temp/*

