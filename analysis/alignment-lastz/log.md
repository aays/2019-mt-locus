log

# 5/11/2018

## making separate fasta files

given mt_loci.fasta:

    from Bio import SeqIO
    d = SeqIO.parse('mt_loci.fasta', 'fasta')
    d = list(d)
    with open('mt_plus.fasta', 'w') as f:
        SeqIO.write(d[0], f, 'fasta')
    with open('mt_minus.fasta', 'w') as f:
        SeqIO.write(d[1], f, 'fasta')

initial lastz run, default parameters: 

    time ./bin/lastz data/references/mt_plus.fasta data/references/mt_minus.fasta \
    --output=data/alignment-lastz/init-out.lav \
    --format=lav

(run took ~49 seconds)

lastz run with different file format:

    time ./bin/lastz data/references/mt_plus.fasta data/references/mt_minus.fasta \
    --output=data/alignment-lastz/init-out.bed \
    --format=general

this format looks much easier to parse in R - we'll be sticking with this

note that it isn't 'officially' bed format - but it's similar enough that I'm calling it that

specification for columns can be found here -
http://www.bx.psu.edu/miller_lab/dist/README.lastz-1.02.00/README.lastz-1.02.00a.html#fmt_general

tomorrow - have a look at the alignment in geneious and eyeball how it lines up


# 7/11/2018

so the lastz flat file output doesn't fully align w/ geneious and I can't quite tell why

geneious alignments seem more contiguous - going to try the gap-free extension method

```bash
time ./bin/lastz data/references/mt_plus.fasta data/references/mt_minus.fasta \
--output=data/alignment-lastz/init-out-gapfree.bed \
--format=general \
--gfextend
```

this resulted in an identical file! seems it's on by default

a gapped alignment?

```bash
time ./bin/lastz data/references/mt_plus.fasta data/references/mt_minus.fasta \
--output=data/alignment-lastz/init-out-gapped.bed \
--format=general \
--gfextend --gapped
```

let's try ungapped in order to avoid the problem of indels

(ldhelmet expects fastas of the same length)

```bash
time ./bin/lastz data/references/mt_plus.fasta data/references/mt_minus.fasta \
--output=data/alignment-lastz/init-out-ungapped.bed \
--format=general --nogapped
```

this looks good - going to stick with it

really small chunks of identically sized matches allows us to avoid that problem almost entirely

manual examination of scores -
- score > 10000 is a good benchmark thus far - but keep manually checking
- few exceptions to this score involve super short tracts
- sometimes tracts overlap slightly but have mismatches at pt of overlap - why?
    - can maybe circumvent this by finding method of joining consecutive high-scoring windows (similar to continuous hotspots)
- method of looking through this in R:

```r
d <- read_tsv('init-out-ungapped.bed')
d %<>% mutate(length1 = end1 - zstart1)
d %<>% tibble::rownames_to_column() # for detecting consecutive windows
d %<>% mutate(idPct = stringr::str_replace(idPct, '%', '')) # convert to numeric
d %<>% mutate(idPct = as.numeric(idPct))
colnames(d)[1] <- 'score' 

ret <- function(position) {
    return(filter(d, zstart1 == position))
}


# 10/11/2018

attempting alignment with 10k score cutoff:

```bash
time ./bin/lastz data/references/mt_plus.fasta data/references/mt_minus.fasta \
--output=data/alignment-lastz/lastz-align-10k-ungapped.bed \
--hspthresh=10000 \
--format=general --nogapped
```

15k?

```bash
time ./bin/lastz data/references/mt_plus.fasta data/references/mt_minus.fasta \
--output=data/alignment-lastz/lastz-align-15k-ungapped.bed \
--hspthresh=15000 \
--format=general --nogapped
```

reciprocal?

```bash
time ./bin/lastz data/references/mt_minus.fasta data/references/mt_plus.fasta \
--output=data/alignment-lastz/lastz-align-10k-ungapped-reciprocal.bed \
--hspthresh=10000 \
--format=general --nogapped
```

next:
- sort through data and find method to deal w/ overlapping windows
- create the relevant fastas w/ vcf2fasta and then query minus sequences based on matches in output table
- _double check the lastz indexing format!_

also, manually checking - it seems that when filtering the previous file for scores > 10000, we have more hits that we do in the newly made file. however, all the hits in the new file are present in the older file, and the ones in the old file that are not in the new one are all of very low scores (see below - using `unique(zstart1)` as a proxy for a unique hit)

```r
d %>% filter(zstart1 %in% d_not) %>% select(score) %>% summary()
     score
 Min.   : 3003
 1st Qu.: 3808
 Median : 4957
 Mean   : 7426
 3rd Qu.:11045
 Max.   :14081
```

# 11/11/2018

## making fastas

from the LD paper - plus strains:

CC-2936
CC-2937
CC-3060
CC-3064
CC-3065
CC-3068
CC-3071
CC-3076
CC-3086

the mt+ locus is from 298298-826737 on chromosome_6

to use Rob's vcf2fasta script, we need to create a symlink to the full quebec vcf

however, in `/scratch/research/projects/chlamydomonas/quebec/all/haplotypeCaller/haplotypeCaller/VCFs` there are now two quebec vcfs - it seems the one without 'all' in the title is missing CC2935, 2936, 2937, and 2938 (ie the Flowers strains)

```python

import vcf
vcf_old = vcf.Reader(filename = 'all_quebec.HC.vcf.gz', compressed = True)
vcf_new = vcf.Reader(filename = 'quebec.HC.vcf.gz', compressed = True)
rec_old = next(vcf_old)
rec_new = next(vcf_new)
samples_old = [c.sample for c in rec_old]
samples_new = [c.sample for c in rec_new]
len(samples_old) # 24
len(samples_new) # 20
set(samples_old).difference(set(samples_new))
# {'CC2936', 'CC2938', 'CC2935', 'CC2937'}

```

after making a symlink to the vcf and vcf.tbi files in data/references, let's make a symlink to vcf2fasta in bin/

```bash
ln -sv /scratch/research/repos/vcf2fasta/vcf2fasta.py bin/
```

the command will be:

```bash
time ./bin/vcf2fasta.py -v data/references/all_quebec.HC.vcf.gz \
-r data/references/mtPlus_ref.chromosome_6.fasta \
-i chromosome_6:298298-826737 \
-s CC2936 CC2937 CC3060 CC3064 CC3065 CC3068 CC3071 CC3076 CC3086 \
--min_GQ 30 > data/aligned-fastas/plus_strains_ref.fasta
```

(took ~4 minutes)

and now for the mt- individuals:

```bash
time ./bin/vcf2fasta.py -v data/references/all_quebec.HC.vcf.gz \
-r data/references/mtMinus_ref.chromosome_6_and_mtMinus.fasta \
-i mtMinus:1-345555 \
-s CC2935 CC2938 CC3059 CC3061 CC3062 CC3063 CC3073 CC3075 CC3079 CC3084 \
--min_GQ 30 > data/aligned-fastas/minus_strains_ref.fasta
```

# 17/11/2018

both the relevant scripts have their first drafts ready:

- analysis/alignment-lastz/clean_lastz_output.R - removes duplicates etc from the lastz output file
- analysis/alignment-lastz/align_mt_fasta.pymt_fasta.py - creates aligned fasta

first pass run:

```bash
Rscript analysis/alignment-lastz/clean_lastz_output.R \
    --file data/alignment-lastz/lastz-align-10k-ungapped.bed \
    --out data/alignment-lastz/lastz-align-filtered.bed
```

and now that this is done, the python script:

```bash
time python3.5 analysis/alignment-lastz/align_mt_fasta.py \
    --plus data/aligned-fastas/plus_strains_ref.fasta \
    --minus data/aligned-fastas/minus_strains_ref.fasta \
    --alignment data/alignment-lastz/lastz-align-filtered.bed \
    --output data/aligned-fastas/mt_aligned.fasta
```

# 18/11/2018

uh oh - we have some overlapping sites.

to fix this, we'll edit the R script, such that:

1. arrange in ascending order of mt+ coordinate
2. if mt+ coordinate at alignment x + 1 is lesser than the ending coordinate of alignment x, this is an overlapping window
3. in these cases, let the end of alignment x be the start of alignment 1 (making x and x + 1 continuous while stil making sure the same sites are being aligned)

so this worked! sort of!

although the fastas seem aligned, the minuses are far longer for some reason:

```python
from Bio import SeqIO
d = [s for s in SeqIO.parse('mt_aligned.fasta', 'fasta')]
[len(str(s.seq)) for s in d]
# minus is 565k, plus is ~520k
```

look into this!

# 19/11/2018

found the issue - the R cleaning script was just making the mt+
coordinates continuous, but not preventing overlap in the mt minus coordinates.

rerunning the commands: 

```bash
Rscript analysis/alignment-lastz/clean_lastz_output.R \
    --file data/alignment-lastz/lastz-align-10k-ungapped.bed \
    --out data/alignment-lastz/lastz-align-filtered.bed
```

the alignment:

```bash
time python3.5 analysis/alignment-lastz/align_mt_fasta.pymt_fasta.py \
    --plus data/aligned-fastas/plus_strains_ref.fasta \
    --minus data/aligned-fastas/minus_strains_ref.fasta \
    --alignment data/alignment-lastz/lastz-align-filtered.bed \
    --output data/aligned-fastas/mt_aligned.fasta
```

# 20/11/2018

although the script is working, it seems upon manual inspection that 
the alignments are all 'shifted over' by one:

```python
>>> for s in plus_refs:
  2     print(s.id, '\n', str(s.seq[500:550]))
CC2936|chromosome_6:298298-826737
 CTGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGTCTGCAGGGGCTGCT
CC2937|chromosome_6:298298-826737
 CTGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGTCTGCAGGGGCTGCT
CC3060|chromosome_6:298298-826737
 CTGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGTTTGCAGGGGCTGCT
CC3064|chromosome_6:298298-826737
 CTGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGTCTGCAGGGGCTGCT
CC3065|chromosome_6:298298-826737
 CTGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGTCTGCAGGGGCTNNN
CC3068|chromosome_6:298298-826737
 CTGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCNNNNTGCNNNNNNNNNN
CC3071|chromosome_6:298298-826737
 CTGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGTCTGCAGGGGCTGCT
CC3076|chromosome_6:298298-826737
 CTGGGCGCGGGCGGGCGGCTGGGCANNNNNNNNNNNNNNNNNNNNNNNNN
CC3086|chromosome_6:298298-826737
 CTGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGTCTGCAGGGGCTGCT
>>> for s in minus_refs:
  2     print(s.id, '\n', str(s.seq[500:550]))
CC2935|mtMinus:1-345555
 NNNNNNNNNNNNNNNNNCTGGGCAGCTCCAGCGGTCTGCAGGGGCTGCTG
CC2938|mtMinus:1-345555
 NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
CC3059|mtMinus:1-345555
 TGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGTCTGCAGGGGCTGCTG
CC3061|mtMinus:1-345555
 TGGGTGCGGGCGGGCGGCTGGGCAGCTCCAGCGGGTTGCAGGGGCTGCTG
CC3062|mtMinus:1-345555
 TGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGTCTGCAGGGGCTGCTG
CC3063|mtMinus:1-345555
 TGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGTCTGCAGGGGCTGCTG
CC3073|mtMinus:1-345555
 TGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGTCTGCAGGGGCTGCTG
CC3075|mtMinus:1-345555
 NNNNNNNNNNNNNNCGGCTGGGCAGCTCCAGCGGGTTGCAGGGGCTGCTG
CC3079|mtMinus:1-345555
 TGGGNGCGGGCGGGCGGCTGGGCAGCTCCAGCGGGTTGCAGGGGCTGCTG
CC3084|mtMinus:1-345555
 TGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGTCTGCAGGGGCTGCTG
>>> show_aln(500, 550)
CC2936|chromosome_6:298298-826737
 CTGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGTCTGCAGGGGCTGCT
CC2937|chromosome_6:298298-826737
 CTGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGTCTGCAGGGGCTGCT
CC3060|chromosome_6:298298-826737
 CTGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGTTTGCAGGGGCTGCT
CC3064|chromosome_6:298298-826737
 CTGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGTCTGCAGGGGCTGCT
CC3065|chromosome_6:298298-826737
 CTGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGTCTGCAGGGGCTNNN
CC3068|chromosome_6:298298-826737
 CTGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCNNNNTGCNNNNNNNNNN
CC3071|chromosome_6:298298-826737
 CTGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGTCTGCAGGGGCTGCT
CC3076|chromosome_6:298298-826737
 CTGGGCGCGGGCGGGCGGCTGGGCANNNNNNNNNNNNNNNNNNNNNNNNN
CC3086|chromosome_6:298298-826737
 CTGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGTCTGCAGGGGCTGCT
CC2938|mtMinus:1-345555
 NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
CC2935|mtMinus:1-345555
 NNNNNNNNNNNNNNNNNCTGGGCAGCTCCAGCGGTCTGCAGGGGCTGCTG
CC3084|mtMinus:1-345555
 TGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGTCTGCAGGGGCTGCTG
CC3059|mtMinus:1-345555
 TGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGTCTGCAGGGGCTGCTG
CC3061|mtMinus:1-345555
 TGGGTGCGGGCGGGCGGCTGGGCAGCTCCAGCGGGTTGCAGGGGCTGCTG
CC3063|mtMinus:1-345555
 TGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGTCTGCAGGGGCTGCTG
CC3073|mtMinus:1-345555
 TGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGTCTGCAGGGGCTGCTG
CC3062|mtMinus:1-345555
 TGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGTCTGCAGGGGCTGCTG
CC3075|mtMinus:1-345555
 NNNNNNNNNNNNNNCGGCTGGGCAGCTCCAGCGGGTTGCAGGGGCTGCTG
CC3079|mtMinus:1-345555
 TGGGNGCGGGCGGGCGGCTGGGCAGCTCCAGCGGGTTGCAGGGGCTGCTG

```

it seems the issue is actually with the plus sequences:

```python
>>> mt_plus_ref = SeqIO.read('../references/mt_plus.fasta', 'fasta')
>>> len(mt_plus_ref)
528439
>>> mt_plus = [s for s in SeqIO.parse('plus_strains_ref.fasta', 'fasta')]
>>> len(mt_plus[0])
528440
>>> mt_plus[0][500:550]
SeqRecord(seq=Seq('CTGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGTCTGCAGGGGCTGCT', SingleLetterAlphabet()), id='CC2936|chromosome_6:298298-826737', name='CC2936|chromosome_6:298298-826737', description='CC2936|chromosome_6:298298-826737', dbxrefs=[])
>>> mt_plus_ref[500:550]
SeqRecord(seq=Seq('TGGGCGCGGGCGGGCGGCTGGGCAGCTCCAGCGGGTTGCAGGGGCTGCTG', SingleLetterAlphabet()), id='chromosome_6', name='chromosome_6', description='chromosome_6', dbxrefs=[])
```

the original vcf2fasta command used python (origin zero) indexing, while vcf2fasta
uses samtools (origin one) indexing. the correct command is therefore

```bash
time ./bin/vcf2fasta.py -v data/references/all_quebec.HC.vcf.gz \
-r data/references/mtPlus_ref.chromosome_6.fasta \
-i chromosome_6:298299-826737 \
-s CC2936 CC2937 CC3060 CC3064 CC3065 CC3068 CC3071 CC3076 CC3086 \
--min_GQ 30 > data/aligned-fastas/plus_strains_ref.fasta
```

third time's the charm?

```bash
time python3.5 analysis/alignment-lastz/align_mt_fasta.pymt_fasta.py \
    --plus data/aligned-fastas/plus_strains_ref.fasta \
    --minus data/aligned-fastas/minus_strains_ref.fasta \
    --alignment data/alignment-lastz/lastz-align-filtered.bed \
    --output data/aligned-fastas/mt_aligned.fasta
```

seems so!

```python
>>> show_aln(2000, 2100)
CC2936|chromosome_6:298299-826737
 AGCAGCAGCTCTAGCAGCACCAGTTGATCAGACAAACGACCGGTCCGTGTGCTGTNNNNNNNNNNNNNNCNCTCCCCGTGTCCTGCTACCTGTTGCGACG
CC2937|chromosome_6:298299-826737
 AGCAGCAGCTCTAGCAGCACCAGTTGATCAGACAAACGACCGGTCCGTGTGCTGTGTGNNNGNNCNNNNCNCTCCCCGTGTCCTGCTACCTGTTGCGACG
CC3060|chromosome_6:298299-826737
 AGCAGCAGCTCTAGCAGCACCAGTTGATCAGACAAACGACCGGTCCGTGTGCTGTGTGACTGACGTTCCCCCTCCCCGTGTCCTGCTACCTGTTGCGACG
CC3064|chromosome_6:298299-826737
 AGCAGCAGCTCTAGCAGCACCAGTTGATCAGACAAACGACCGGTCCGTGTGCTGTGTGACTGACGTTCCCCCTCCCCGTGTCCTGCTACCTGTTGCGACG
CC3065|chromosome_6:298299-826737
 NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
CC3068|chromosome_6:298299-826737
 AGCAGCAGCTCTAGCAGCACCAGTTGATCAGACAAACGACCGGTCCGTGTGCTGTGNNNNNNNNNNNNNCNCTCCCCGTGTCCTGCTACCTGTTGCGACG
CC3071|chromosome_6:298299-826737
 NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
CC3076|chromosome_6:298299-826737
 NNCAGCAGCTCTAGCAGCACCAGTTGATCAGACAAACGACCGGTCCGTGNNNNNNNNNNNNNNNNNNNNCNCTCCCCGTGTCCTGCTACCTGTTGCGACG
CC3086|chromosome_6:298299-826737
 AGCAGCAGCTCTAGCAGCACCAGTTGATCAGACAAACGACCGGTCCGTGTGCTGTGNNNNNNNNNNNNNCNCTCCCCGTGTCCTGCTACCTGTTGCGACG
CC3075|mtMinus:1-345555
 NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
CC2935|mtMinus:1-345555
 AGCAGCAGCTCTAGCAGCACCAGTTGATCAGANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTCCCCNTNCCCGTGTCCTGCTACCTGTTGCGACG
CC3073|mtMinus:1-345555
 NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCNNNNNNNNTGTCCTGCTACCTGTTGCGACG
CC3084|mtMinus:1-345555
 AGCAGCAGCTCTAGCAGCACCAGTTGATCAGACAAACGACCGGTCCGTGTGCTGTGTGANNNNNNNNNNCNCTCCCCGTGTCCTGCTACCTGTTGCGACG
CC3059|mtMinus:1-345555
 AGCAGCAGCTCTAGCAGCACCAGTTGATCAGACAAACGACCGGTCCGTGTGCTGTGTGANNNNNNNNNNCNCTCCCCGTGTCCTGCTACCTGTTGCGACG
CC3079|mtMinus:1-345555
 NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
CC3063|mtMinus:1-345555
 AGCAGCAGCTCTAGCAGCACCAGTTGATCAGACAAACGACCGGTCCGTGTGCTGTGTGANNNNNNNNNNCNCTCCCCGTGTCCTGCTACCTGTTGCGACG
CC3062|mtMinus:1-345555
 NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCNNNNCCCGTGTCCTGCTACCTGTTGCGACG
CC3061|mtMinus:1-345555
 AGCAGCAGCTGTAGCAGCACCAGTTGATCAGACAAACGACCGGTCCGTGTGCTGTGTGACTGACCTTCCCCTTACCCGTGTCCTGCTACCTGTTGCGACG
CC2938|mtMinus:1-345555
 AGCAGCAGCTCTAGCAGCACCAGTTGATCAGACAAACGACCGGTCCGTGTGCTGTGTGANNNNNNNNNNCNNNNNCCGTGTCCTGCTACCTGTTGCGACG
```

onto using LDhelmet to calculate recombination rates!












