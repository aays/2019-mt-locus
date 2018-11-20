mkdir data/recombination-ldhelmet/intermediate-files
mkdir data/recombination-ldhelmet/recombination-estimates

time ./bin/ldhelmet find_confs \
--num_threads 10 \
--window_size 50 \
--output_file data/recombination-ldhelmet/intermediate-files/output.conf \
data/aligned-fastas/mt_aligned.fasta

time ./bin/ldhelmet table_gen \
--num_threads 10 \
--conf_file data/recombination-ldhelmet/intermediate-files/output.conf \
--theta 0.01 \
--rhos 0.0 0.1 10.0 1.0 100.0 \
--output_file data/recombination-ldhelmet/intermediate-files/output.lk

time ./bin/ldhelmet pade \
--num_threads 10 \
--conf_file data/recombination-ldhelmet/intermediate-files/output.conf \
--theta 0.01 \
--output_file data/recombination-ldhelmet/intermediate-files/output.pade

time ./bin/ldhelmet rjmcmc \
--num_threads 10 \
--window_size 50 \
--seq_file data/aligned-fastas/mt_aligned.fasta \
--lk_file data/recombination-ldhelmet/intermediate-files/output.lk \
--pade_file data/recombination-ldhelmet/intermediate-files/output.pade \
--num_iter 1000000 \
--burn_in 100000 \
--block_penalty 50 \
--output_file data/recombination-ldhelmet/intermediate-files/output.post

time ./bin/ldhelmet post_to_text \
--mean \
--perc 0.025 \
--perc 0.50 \
--perc 0.975 \
--output_file data/recombination-ldhelmet/recombination-estimates/mt_locus_recombination.txt
data/recombination-ldhelmet/intermediate-files/output.post
