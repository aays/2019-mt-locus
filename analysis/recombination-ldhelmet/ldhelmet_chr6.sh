time ./bin/ldhelmet find_confs \
--num_threads 10 \
--window_size 50 \
--output_file data/recombination-ldhelmet/intermediate-files/chr6.conf \
data/aligned-fastas/chromosome_6_all.fasta

time ./bin/ldhelmet table_gen \
--num_threads 10 \
--conf_file data/recombination-ldhelmet/intermediate-files/chr6.conf \
--theta 0.01 \
--rhos 0.0 0.1 10.0 1.0 100.0 \
--output_file data/recombination-ldhelmet/intermediate-files/chr6.lk

time ./bin/ldhelmet pade \
--num_threads 10 \
--conf_file data/recombination-ldhelmet/intermediate-files/chr6.conf \
--theta 0.01 \
--output_file data/recombination-ldhelmet/intermediate-files/chr6.pade

time ./bin/ldhelmet rjmcmc \
--num_threads 10 \
--window_size 50 \
--seq_file data/aligned-fastas/chromosome_6_all.fasta \
--lk_file data/recombination-ldhelmet/intermediate-files/chr6.lk \
--pade_file data/recombination-ldhelmet/intermediate-files/chr6.pade \
--num_iter 1000000 \
--burn_in 100000 \
--block_penalty 50 \
--output_file data/recombination-ldhelmet/intermediate-files/chr6.post

time ./bin/ldhelmet post_to_text \
--mean \
--perc 0.025 \
--perc 0.50 \
--perc 0.975 \
--output_file data/recombination-ldhelmet/recombination-estimates/chr6_all_recombination.txt \
data/recombination-ldhelmet/intermediate-files/chr6.post
