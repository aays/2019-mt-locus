# run from root dir!

mkdir -p data/recombination-ldhelmet/intermediate-files
mkdir -p data/recombination-ldhelmet/recombination-estimates

count=1
filecount="$(ls data/aligned-fastas/alignments/ | wc -l)"

for filename in data/aligned-fastas/alignments/*fasta; do
    base=$(basename $filename .fasta)

    echo "LDhelmet run for ${base}"
    echo "File ${count} of ${filecount}"

    time ./bin/ldhelmet find_confs \
    --num_threads 10 \
    --window_size 50 \
    --output_file data/recombination-ldhelmet/intermediate-files/${base}.conf ${infile}

    time ./bin/ldhelmet table_gen \
    --num_threads 10 \
    --conf_file data/recombination-ldhelmet/intermediate-files/${base}.conf \
    --theta 0.01 \
    --rhos 0.0 0.1 10.0 1.0 100.0 \
    --output_file data/recombination-ldhelmet/intermediate-files/${base}.lk

    time ./bin/ldhelmet pade \
    --num_threads 10 \
    --conf_file data/recombination-ldhelmet/intermediate-files/${base}.conf \
    --theta 0.01 \
    --output_file data/recombination-ldhelmet/intermediate-files/${base}.pade

    time ./bin/ldhelmet rjmcmc \
    --num_threads 10 \
    --window_size 50 \
    --seq_file ${filename} \
    --lk_file data/recombination-ldhelmet/intermediate-files/${base}.lk \
    --pade_file data/recombination-ldhelmet/intermediate-files/${base}.pade \
    --num_iter 1000000 \
    --burn_in 100000 \
    --block_penalty 50 \
    --output_file data/recombination-ldhelmet/intermediate-files/${base}.post

    time ./bin/ldhelmet post_to_text \
    --mean \
    --perc 0.025 \
    --perc 0.50 \
    --perc 0.975 \
    --output_file data/recombination-ldhelmet/recombination-estimates/${base}.txt \
    data/recombination-ldhelmet/intermediate-files/${base}.post

    # increment counter
    (( count ++ ))

    echo "Removing temp files..."
    rm -v data/recombination-ldhelmet/intermediate-files/${base}.*
