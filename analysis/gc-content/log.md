
## 18/1/2019

overall to do: calculate GC content at 4D sites in non-gametologous regions

things needed:
- plus NG fasta (in `aligned-fastas`)
- annotation table parser - 4D sites defined as intronic/intergenic/4D annotations in there
- script to transpose fasta? the coordinates should match the mt+ though... might not need transposition here
- script that iterates through alignment one base at a time, checks 4D status, and calculates GC if true

## 22/1/2019

wait, are we calculating GC from the reference? or from the consensus of the plus sequences?

in either case, let's get the annotation table infrastructure going:

```bash
cd data/references
ln -sv /scratch/research/references/chlamydomonas/5.3_chlamy_w_organelles_mt_minus/annotation/concatenated_GFF/annotation_table.txt.gz* .
cd ../../
cp -v ../genomewide_recombination/analysis/fiftyshadesofgreen/annotation_parser/ant.py analysis/gc-content/
```

let's start the script with the assumption we're calculating GC content from the consensus.
this script can iterate through the NG fasta given an offset value, ignoring any sites that are Ns,
and returning GC content in a given windowsize. other things to return would be

1. total usable sites
2. count of Gs and Cs

so that overall GC content in larger windows (or even across all NG regions)
can be recalculated from this outfile.

testing the script:

```bash
time python3.5 analysis/gc-content/gc_calc.py \
--filename data/aligned-fastas/plus_non_gametolog.fasta \
--annotation data/references/annotation_table.txt.gz \
--windowsize 1000 \
--region chromosome_6:298299-826737 \
--outfile data/gc-content/gc_4D_NG_1k.txt
```

seems this GC content level is just below genomewide levels:

```R
> d %>% summarise(total_GC = sum(GC), total_sites = sum(total_sites)) %>% 
+ mutate(GC_content = total_GC / total_sites)
# A tibble: 1 x 3
  total_GC total_sites GC_content
     <int>       <int>      <dbl>
1    74828      123494  0.6059242
```

