# ldhelmet on the non-gametolog parts of the mt alleles individually

for allele in plus minus; do
    filename=${allele}_non_gametolog.fasta;

    time ./bin/ldhelmet find_confs \
    --num_threads 10 \
    --window_size 50 \
    --output_file data/recombination-ldhelmet/intermediate-files/${allele}.conf \
    data/aligned-fastas/${filename}

    time ./bin/ldhelmet table_gen \
    --num_threads 10 \
    --conf_file data/recombination-ldhelmet/intermediate-files/${allele}.conf \
    --theta 0.01 \
    --rhos 0.0 0.1 10.0 1.0 100.0 \
    --output_file data/recombination-ldhelmet/intermediate-files/${allele}.lk

    time ./bin/ldhelmet pade \
    --num_threads 10 \
    --conf_file data/recombination-ldhelmet/intermediate-files/${allele}.conf \
    --theta 0.01 \
    --output_file data/recombination-ldhelmet/intermediate-files/${allele}.pade

    time ./bin/ldhelmet rjmcmc \
    --num_threads 10 \
    --window_size 50 \
    --seq_file data/aligned-fastas/${filename} \
    --lk_file data/recombination-ldhelmet/intermediate-files/${allele}.lk \
    --pade_file data/recombination-ldhelmet/intermediate-files/${allele}.pade \
    --num_iter 1000000 \
    --burn_in 100000 \
    --block_penalty 50 \
    --output_file data/recombination-ldhelmet/intermediate-files/${allele}.post

    time ./bin/ldhelmet post_to_text \
    --mean \
    --perc 0.025 \
    --perc 0.50 \
    --perc 0.975 \
    --output_file data/recombination-ldhelmet/recombination-estimates/${allele}_only_recombination.txt \
    data/recombination-ldhelmet/intermediate-files/${allele}.post;

done
