log for alignment analysis

## 5/11/2018

### making separate fasta files

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


## 7/11/2018

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
```

## 10/11/2018

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

## 11/11/2018

### making fastas

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

## 17/11/2018

both the relevant scripts have their first drafts ready:

- `analysis/alignment-lastz/clean_lastz_output.R` - removes duplicates etc from the lastz output file
- `analysis/alignment-lastz/align_mt_fasta.py` - creates aligned fasta

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

## 18/11/2018

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

## 19/11/2018

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

## 20/11/2018

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
time python3.5 analysis/alignment-lastz/align_mt_fasta.py \
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

## 24/11/2018

I just realized a significant bug - the fact that the aligned fasta file contains all of the
mt+ sequence but large chunks of missing mt- sequence means that the LDhelmet output there is 
not a 'true' representation of RR in gametologs - LDhelmet is still calculating 
mt+ only 'recombination' in those non-shared regions. the script should instead create a
fasta that has _only_ gametolog regions.

now that the script has been edited, let's try this again:

```bash
time python3.5 analysis/alignment-lastz/align_mt_fasta.py \
    --plus data/aligned-fastas/plus_strains_ref.fasta \
    --minus data/aligned-fastas/minus_strains_ref.fasta \
    --alignment data/alignment-lastz/lastz-align-filtered.bed \
    --output data/aligned-fastas/mt_aligned_only.fasta
```

seems to have worked!

there's a gap in ~1415-1450 - is that reflected in the fasta?

```python
>>> from Bio import SeqIO
>>> d = [s for s in SeqIO.parse('mt_aligned_only.fasta', 'fasta')]
>>> def show_aln(start, end):
        for s in d:
            print(s.id, '\n', s.seq)
>>> show_aln(1410, 1460)
CC3064|chromosome_6:298299-826737
 TTGCATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCCGTGACCCCC
CC3071|chromosome_6:298299-826737
 TTGCATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCCGTGACCCCC
CC3076|chromosome_6:298299-826737
 TTGCATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCCGTGACCCCC
CC2936|chromosome_6:298299-826737
 TTGCATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCCGTGACCCCC
CC3068|chromosome_6:298299-826737
 TTGCATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCCGTGACCCCC
CC3060|chromosome_6:298299-826737
 TTGCATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCCGTGACCCCC
CC3065|chromosome_6:298299-826737
 TTGCATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCCGTGACCCCC
CC2937|chromosome_6:298299-826737
 TTGCATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCCGTGACCCCC
CC3086|chromosome_6:298299-826737
 TTGCATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCCGTGACCCCC
CC3073|mtMinus:1-345555
 NNNNANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
CC2935|mtMinus:1-345555
 NNNNANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCCC
CC3063|mtMinus:1-345555
 NNNNANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
CC3084|mtMinus:1-345555
 NNNNANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
CC3079|mtMinus:1-345555
 TTGCATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCCGTTACCCCC
CC3062|mtMinus:1-345555
 NNNNANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
CC3061|mtMinus:1-345555
 TTGCATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCCGTGACCCCC
CC3059|mtMinus:1-345555
 NNNNANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
CC3075|mtMinus:1-345555
 TTGCATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCCGTTACCCCC
CC2938|mtMinus:1-345555
 NNNNANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCCGTGACCCCC

```

## 25/11/2018

so apparently I've been using the wrong VCFs - hence the large amounts of
missing sequence in the alignment above.

and so we vcf2fasta again:

```bash
time ./bin/vcf2fasta.py -v data/references/all_quebec.mtPlus.gGVCFs.vcf.gz \
-r data/references/mtPlus_ref.chromosome_6.fasta \
-i chromosome_6:298299-826737 \
-s CC2936 CC2937 CC3060 CC3064 CC3065 CC3068 CC3071 CC3076 CC3086 \
--min_GQ 30 > data/aligned-fastas/plus_strains_ref.fasta

time ./bin/vcf2fasta.py -v data/references/all_quebec.mtMinus.gGVCFs.vcf.gz \
-r data/references/mtMinus_ref.chromosome_6_and_mtMinus.fasta \
-i mtMinus:1-345555 \
-s CC2935 CC2938 CC3059 CC3061 CC3062 CC3063 CC3073 CC3075 CC3079 CC3084 \
--min_GQ 30 > data/aligned-fastas/minus_strains_ref.fasta
```

and now for the alignment (again...)

```bash
time python3.5 analysis/alignment-lastz/align_mt_fasta.py \
    --plus data/aligned-fastas/plus_strains_ref.fasta \
    --minus data/aligned-fastas/minus_strains_ref.fasta \
    --alignment data/alignment-lastz/lastz-align-filtered.bed \
    --output data/aligned-fastas/mt_aligned_only.fasta
```

## 28/11/2018

things to do after meeting with Rob:
- what's the difference in lastz output between gapped and nongapped alignments?
    - ie how much data, if any, are we losing with the gap free method?
- make fastas for the plus and minus alleles that specifically do _not_ contain gametologs
- write a script that correctly calculates the weighted means, and divides by the correct amt of sequence
    - the current R script divides by the entire length of the mt locus, which is incorrect

### lastz - gapped versus ungapped

we can start by running the exact same lastz command as before, minus the gapfree arg

```bash
time ./bin/lastz data/references/mt_plus.fasta data/references/mt_minus.fasta \
--output=data/alignment-lastz/lastz-align-10k-gapped.bed \
--hspthresh=10000 \
--format=general
```

off the bat, it seems we get less hits, but that the hits are much larger in the gapped version:

```bash
$ wc -l *
   261 lastz-align-10k-gapped.bed
  1004 lastz-align-10k-ungapped.bed
```

## 29/11/2018

so whether or not we use gaps (that analysis is running right now and being logged in 
`analysis/recombination-ldhelmet/log.md`), we'll be needing a script that generates fastas
containing non-gametolog regions.

```bash
cp align_mt_fasta.py make_mt_only.py
```

### update

now that this script has been written, let's give it a shot
and see if there are any outstanding errors:

```bash
Rscript analysis/alignment-lastz/clean_lastz_output.R \
    --file data/alignment-lastz/lastz-align-10k-gapped.bed \
    --outfile data/alignment-lastz/lastz-align-gapped-filtered.bed

time python3.5 analysis/alignment-lastz/make_mt_only.py \
    --fasta data/aligned-fastas/minus_strains_ref.fasta \
    --alignment data/alignment-lastz/lastz-align-gapped-filtered.bed \
    --mt_allele minus \
    --output data/alignment-lastz/minus_gametolog_masked.fasta
```

okay so this may have not gone according to plan

```python
>>> for s in SeqIO.parse(filename, 'fasta'):
  2     print(str(s.seq).count('N'), len(str(s.seq)))
345555 345555
345555 345555
345555 345555
345555 345555
345555 345555
345555 345555
345555 345555
345555 345555
345555 345555
345555 345555
```

first, do we see the same thing happening w/ the plus?
(ie does this have to do with the rev comp business?)

```bash
time python3.5 analysis/alignment-lastz/make_mt_only.py \
    --fasta data/aligned-fastas/minus_strains_ref.fasta \
    --alignment data/alignment-lastz/lastz-align-gapped-filtered.bed \
    --mt_allele minus \
    --output data/alignment-lastz/minus_gametolog_masked.fasta
```

```python
>>> from Bio import SeqIO
>>> filename = 'plus_gametolog_masked.fasta'
>>> for s in SeqIO.parse(filename, 'fasta'):
  2     print(s.seq.count('N'), len(s.seq))
528439 528439
528439 528439
528439 528439
528439 528439
528439 528439
528439 528439
528439 528439
528439 528439
528439 528439
```

look into fixing this tomorrow.

in the meantime, let's give the original alignment script a shot using the
gapped alignment:

```bash
time python3.5 analysis/alignment-lastz/align_mt_fasta.py \
    --plus data/aligned-fastas/plus_strains_ref.fasta \
    --minus data/aligned-fastas/minus_strains_ref.fasta \
    --alignment data/alignment-lastz/lastz-align-gapped-filtered.bed \
    --output data/aligned-fastas/mt_aligned_gapped.fasta
```

while this didn't raise any errors, it seems the
minus sequence is ~100 kb shorter than the plus:

```python
>>> for s in SeqIO.parse(filename, 'fasta'):
  2     print(s.id, '\n', len(s.seq))
CC3060|chromosome_6:298299-826737
 528439
CC3086|chromosome_6:298299-826737
 528439
CC3064|chromosome_6:298299-826737
 528439
CC3065|chromosome_6:298299-826737
 528439
CC3071|chromosome_6:298299-826737
 528439
CC2936|chromosome_6:298299-826737
 528439
CC2937|chromosome_6:298299-826737
 528439
CC3076|chromosome_6:298299-826737
 528439
CC3068|chromosome_6:298299-826737
 528439
CC3079|mtMinus:1-345555
 528396
CC3061|mtMinus:1-345555
 528396
CC3073|mtMinus:1-345555
 528396
CC2935|mtMinus:1-345555
 528396
CC3062|mtMinus:1-345555
 528396
CC3075|mtMinus:1-345555
 528396
CC3063|mtMinus:1-345555
 528396
CC3059|mtMinus:1-345555
 528396
CC2938|mtMinus:1-345555
 528396
CC3084|mtMinus:1-345555
 528396
```

looking over at the very end, they don't seem especially aligned:

```python
>>> for s in SeqIO.parse(filename, 'fasta'):
  2     print(s.id, '\n', str(s.seq)[528300:528396])
CC3060|chromosome_6:298299-826737
 CTTGTTAACGGTATGACGCCACACGAGCGAGACTGTAGGCGCCATTTGCAGCTTTATCTGACAAATACCGACAGAAGACNTGCGGCGACTGAAGTT
CC3086|chromosome_6:298299-826737
 CTTGTTAACGGTATGACGCCACACGAGCGAGACTGTAGGCGCCATTTGCAGCTGTATCTGACAAATACCGACAGAAGACTTGCGGCGACTGAAGTT
CC3064|chromosome_6:298299-826737
 CTTGTTAACGGTATGACGCCACACGAGCGAGACTGTAGGCGCCATTTGCAGCTTTATCTGACAAATACCGACAGAAGACNTGCGGCGACTGAAGTT
CC3065|chromosome_6:298299-826737
 CTTGTTAACGGTATGACGCCACACGAGCGAGACTGTAGGCGCCATTTGCAGCTTTATCTGACAAATACCGACAGAAGACNTGCGGCGACTGAAGTT
CC3071|chromosome_6:298299-826737
 CTTGTTAACGGTATGACGCCACACGAGCGAGACTGTAGGCGCCATTTGCAGCTTTATCTGACAAATACCGACAGAAGACNTGCGGCGACTGAAGTT
CC2936|chromosome_6:298299-826737
 CTTGTTAACGGTATGACGCCACACGAGCGAGACTGTAGGCGCCATTTGCAGCTTTATCTGACAAATACCGACAGAAGACTTGCGGCGACTGAAGTT
CC2937|chromosome_6:298299-826737
 CTTGTTAACGGTATGACGCCACACGAGCGAGACTGTAGGCGCCATTTGCAGCTTTATCTGACAAATACCGACAGAAGACTTGCGGCGACTGAAGTT
CC3076|chromosome_6:298299-826737
 CTTGTTAACGGTATGACGCCACACGAGCGAGACTGTAGGCGCCATTTGCAGCTTTATCTGACAAATACCGACNGNANNCTTGCGGCGACTGAAGTT
CC3068|chromosome_6:298299-826737
 CTTGTTAACGGTATGACGCCACACGAGCGAGACTGTAGGCGCCATTTGCAGCTTTATCTGACAAATACCGACAGAAGACNNNNNNNGACTGAAGTT
CC3079|mtMinus:1-345555
 NNNNNNNNGNNNNCTATGCACGCACCAGGCGTCCATGTTGCGGCGACTGAAGTTGATNNNNNNGNNNNNGNNNNNNGNNNNNNNNNNNNNNNNNNN
CC3061|mtMinus:1-345555
 CANNNNNNGNNNNCTATGGACGCACCAGGCGTCCATGTTGCGGCGACTGAAGTTGATNNNNNNGNNNNNGNNCNCNGNNNNNNNNNNNNNNNNNNN
CC3073|mtMinus:1-345555
 CAGNNNNNGNNNNCTATGCACGCACCAGGCGTCCATGTTGCGGCGACTGAAGTTGATNNNNNNGNNNNNGNNNNNNGNNNNNNNNNNNNNNNNNNN
CC2935|mtMinus:1-345555
 CAGNNNACGNNNNCTATGCACGCACCAGGCGTCCATGTTGCGGCGACTGAAGTTGATNNNNNNGNNNNNGNNNNNNGNNNNNNNNNNNNNNNNNNN
CC3062|mtMinus:1-345555
 CAGNNNNNGNNNNCTATGCACGCACCAGGCGTCCATGTTGCGGCGACTGAAGTTGATNNNNNNGNNNNNGNNNNNNGNNNNNNNNNNNNNNNNNNN
CC3075|mtMinus:1-345555
 CAGAAGACGNNNNCTATGCACGCACCAGGCGTCCATGTTGCGGCGACTGAAGTTGATNNNNNNGNNNNNGNNNNNNGNNNNNNNNNNNNNNNNNNN
CC3063|mtMinus:1-345555
 CAGAAGACGNNNNCTATGCACGCACCAGGCGTCCATGTTGCGGCGACTGAAGTTGATNNNNNNGNNNNNGNNNNNNGNNNNNNNNNNNNNNNNNNN
CC3059|mtMinus:1-345555
 CAGNNNNNGNNNNCTATGCACGCACCAGGCGTCCATGTTGCGGCGACTGAAGTTGATNNNNNNGNNNNNGNNNNNNGNNNNNNNNNNNNNNNNNNN
CC2938|mtMinus:1-345555
 CAGAAGACGACNNCTATGCACGCACCAGGCGTCCATGTTGCGGCGACTGAAGTTGATNNNNNNGNNNNNGNNNNNNGNNNNNNNNNNNNNNNNNNN
CC3084|mtMinus:1-345555
 CAGNNGACGNNNNCTATGCACGCACCAGGCGTCCATGTTGCGGCGACTGAAGTTGATNNNNNNGNNNNNGNNNNNNGNNNNNNNNNNNNNNNNNNN
```

I think this is because alignments aren't the same length:

```R
> d <- read_tsv('lastz-align-gapped-filtered.bed', col_types = cols())
> d %>% mutate(length1 = end1 - zstart1, length2 = end2 - zstart2) %>%
+ filter(length1 != length2) %>%
+ dim()
[1]  0 17
```

guess not.

another thing to look into tomorrow!

in the meantime, I've noticed that the scores are MUCH higher
in the gapped alignment.

```R
> colnames(d)[1] <- 'score'
> summary(d$score)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  10709   20791   37506  147630   89860 3857474
```

as compared to:

```R
> colnames(d_old)[1] <- 'score'
> summary(d_old$score)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  10024   13995   20237   32111   35346  290108
```

what does a higher threshold give us?

```bash
time ./bin/lastz data/references/mt_plus.fasta data/references/mt_minus.fasta \
--output=data/alignment-lastz/lastz-align-100k-gapped.bed \
--hspthresh=100000 \
--format=general
```

okay, doesn't really pass the eye test - clear gametologs from the 
geneious alignment have been excluded. 

## 30/11/2018

today:
- debugging the allele separated script - why are we getting all Ns?
- debugging the main alignment script - why are the two strains diff lengths?

### allele separated script

if I have the script return the number of non shared bases
vs the total length of the sequence, I get 245406 unique bases
and the expected seq length of 528409 - which tells me something's
going wrong in the output sequence generation

update: the script was masking stuff but then just rewriting the input file
back to the outfile! unbelievable.

```bash
time python3.5 analysis/alignment-lastz/make_mt_only.py \
--fasta data/aligned-fastas/plus_strains_ref.fasta \
--alignment data/alignment-lastz/lastz-align-gapped-filtered.bed \
--mt_allele plus \
--output data/aligned-fastas/plus_non_gametolog.fasta

time python3.5 analysis/alignment-lastz/make_mt_only.py \
--fasta data/aligned-fastas/minus_strains_ref.fasta \
--alignment data/alignment-lastz/lastz-align-gapped-filtered.bed \
--mt_allele minus \
--output data/aligned-fastas/minus_non_gametolog.fasta
```

### main aligned script

so the current issue is that strains from the two mating types seem
to be slightly misaligned, and that the plus side is longer than the minus

first, are there size discrepancies in the filtered gapped alignment itself?

```R
> d <- read_tsv('lastz-align-gapped-filtered.bed', col_types = cols())
> colnames(d)[1] <- 'score'
> d %>%
+ mutate(length1 = end1 - zstart1, length2 = end2 - zstart2) %>%
+ summarise(plus_sum = sum(length1), minus_sum = sum(length2))
# A tibble: 1 x 2
  plus_sum minus_sum
     <int>     <int>
1   283033    283033
```

doesn't seem to be the case - but despite this, the minus alignment is
shorter than the plus

```bash
time python3.5 analysis/alignment-lastz/align_mt_fasta.py \
--plus data/aligned-fastas/plus_strains_ref.fasta \
--minus data/aligned-fastas/minus_strains_ref.fasta \
--alignment data/alignment-lastz/lastz-align-gapped-filtered.bed \
--output test.out
```

## 1/12/2018

testing alternate lastz output formats to
examine the gaps:

```bash
time ./bin/lastz data/references/mt_plus.fasta data/references/mt_minus.fasta \
--output=data/alignment-lastz/lastz-align-10k-gapped.maf \
--hspthresh=10000 \
--format=maf
```

let's test whether this aligns with the nogap thing:

```bash
sed -n '16,17p' lastz-align-10k-gapped.maf > test_align.maf
```

not quite, upon manual inspection. we'll have to iterate through
the maf file and create a new mt aligning script

## 3/12/2018

so this is proving a bit difficult - I'm updating the script to read the bed as
an ancillary file, which it will use to avoid duplicates

cleaning the bed:

```R
library(tidyverse)
library(magrittr)
d <- read_tsv('lastz-align-10k-gapped.bed')
colnames(d)[1] <- 'score'
d %<>%
    group_by(zstart1) %>%
    filter(score == max(score)) %>% 
    filter(zstart2 == min(zstart2)) %>%
    ungroup()
write_tsv(d, 'lastz-align-10k-gapped-filtered.bed')
```

instead of the original single-fasta plan:

- from each maf alignment, if matched in filtered bed, create a fasta file containing the alignment
    - filtered bed will have removed multiple mt- matches and picked the one with the highest score
    - this alignment will grab the relevant mt+ and mt- strains
    - the header of the fasta should contain the _actual_ coordinates of the sequence
        - by actual, this should mean both the coordinates wrt the mt AND chromosome 6 as a whole
    - due to the gapped alignment these will likely not correspond to the length of the seqs
- next, run LDhelmet on each of these to get recombination rate estimates for each alignment
- creating a 'full mt locus fasta':
    - this will follow the mt+ coordinates, as originally planned
    - however, handling of gaps will be different
        - in cases where mt+ has a gap, the mt- sequence at the gap will be omitted
        - ie insertions in the mt- will be removed
        - if the mt- has a gap, we just ignore it/replace it with Ns

note that in the maf file, the size reported does _not_ contain gaps,
nor does the interval size in the equivalent bed file.


## 12/12/2018

first pass attempt at the script:

```bash
mkdir data/aligned-fastas/alignments

time python3.5 analysis/alignment-lastz/align_mt_fasta_maf.py \
--plus data/aligned-fastas/plus_strains_ref.fasta \
--minus data/aligned-fastas/minus_strains_ref.fasta \
--alignment data/alignment-lastz/lastz-align-10k-gapped.maf \
--bed data/alignment-lastz/lastz-align-10k-gapped-filtered.bed \
--outdir data/aligned-fastas/alignments
```

it works! and it only takes 9 seconds!!!

next - scope out some of these manually in the maf
to make sure they were correctly added - ie the one at 760104


## 13/12/2018

upon manual inspection, it's looking like some of these were
actually misaligned, and that the same chunk of strain
sequence keeps getting drawn from the references

the issue looks to be that when the script draws from
the strain sequences, it always starts from 0, and not
where the actual match starts

trying the script again after fixing this:

```
time python3.5 analysis/alignment-lastz/align_mt_fasta_maf.py \
--plus data/aligned-fastas/plus_strains_ref.fasta \
--minus data/aligned-fastas/minus_strains_ref.fasta \
--alignment data/alignment-lastz/lastz-align-10k-gapped.maf \
--bed data/alignment-lastz/lastz-align-10k-gapped-filtered.bed \
--outdir data/aligned-fastas/alignments/
```

this now raises an `IndexError` - it seems there is
an `i + start1` combination that's somehow larger than
the entire sequence - check the ones at the end
2. when looking at the ones at the end, make sure the
maf actually aligns with what's in the bed - right now,
it sort of looks like there are some extra files in the output
(although they're listed earlier in the bed... why are these
out of order again? maybe remake the bed so that it's in
ascending order of `zstart1` for easier comparison)

## 14/12/2018

THE SCRIPT WORKED!

although I'm noticing from the output files that
some have different starts but end at the same position (ie
ending at 796207)

though that being said, LDhelmet will likely calculate 
the same values for the same regions in both - and so
these can be concatenated in the LDhelmet outputs
and duplicates removed

finally, let's make the mt-separated versions of these
strains using the new filtered bed file:

```bash
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
```

for some reason this didn't work properly on the minus (based on the eye test)

## 16/12/2018

today: debugging what went wrong with the minus non-gametolog file

update - not sure what was going wrong - manually reran the functions in the
script and everything added up, and upon rerunning the script the correct
output file was generated. must have made a typo the first time...

onto LDhelmet for these files!

however, it's worth keeping in mind that later on, we'll
be needing to create a single mt-aligned file that follows the mt+
coordinates as a reference

this file will basically ignore gaps in the mt+ and only keep
gaps in the mt-, so as to preserve the mt+ coordinates

## 18/12/2018

new alignment method (motivation of problem discussed in `recombination-ldhelmet` log)

so:

`transpose_fastas.py` -
1. open transposed file to write to
2. for each fasta in `alignments`:
    a. for each position in the fasta:
        a. if mt+ is not '-': write out position and base for strains
        b. elif mt+ is '-': skip over this site

although to make sure we get the strain names,
we have to modify `align_mt_fasta_maf.py` to output
those as fasta ids (instead of just the regions)

```
mkdir -p data/aligned-fastas/alignments_strains
time python3.5 analysis/alignment-lastz/align_mt_fasta_maf.py \
--plus data/aligned-fastas/plus_strains_ref.fasta \
--minus data/aligned-fastas/minus_strains_ref.fasta \
--alignment data/alignment-lastz/lastz-align-10k-gapped.maf \
--bed data/alignment-lastz/lastz-align-10k-gapped-filtered.bed \
--outdir data/aligned-fastas/alignments_strains

time python3.5 analysis/alignment-lastz/transpose_fastas.py \
--directory data/aligned-fastas/alignments_strains/ \
--outfile data/aligned-fastas/mt_aligned_transposed.txt
```

`remove_duplicates.R` 
1. read in file above as df
2. dplyr::distinct() to remove duplicates
3. assert that no position appears twice

```bash
Rscript analysis/alignment-lastz/remove_duplicates.R \
data/aligned-fastas/mt_aligned_transposed.txt \
data/aligned-fastas/mt_aligned_transposed_filtered.txt
```

and then:

`combine_fastas.py` - this input file will have no mt+ gaps or duplicates
1. open transposed file
2. create dict with strain names from header
3. instantiate counter = 298299
3. for line in fasta:
    a. check that position = counter
    b. grow sequence string (dict value) for all strains
    c. if not position = counter # gap happened
        a. gap_size = position - counter
        b. increment sequence string with gap_size number of Ns
        c. then proceed to next iteration

```bash
time python3.5 analysis/alignment-lastz/combine_fastas.py \
--file data/aligned-fastas/mt_aligned_transposed_filtered.txt \
--outfile data/aligned-fastas/mt_aligned_final.fasta
```


## 20/12/2018

lastz alignment with mt- as the target -

here, we don't necessarily need the gaps - the bed coordinates
are sufficient to denote regions that are (and aren't) gametologous

after running this, we need to
1. check whether the 'coverage' of aligned bases across the two mt alleles has changes
2. if so, then updating scripts that mark non-gametologous regions with both sets

lastz run:

```bash
time ./bin/lastz data/references/mt_minus.fasta data/references/mt_plus.fasta \
--output=data/alignment-lastz/lastz-align-10k-gapped-reciprocal.bed \
--hspthresh=10000 \
--format=general
```

filtering:

```R
library(tidyverse)
library(magrittr)
d <- read_tsv('lastz-align-10k-gapped-reciprocal.bed')
colnames(d)[1] <- 'score'
d %<>%
    group_by(zstart1) %>%
    filter(score == max(score)) %>% 
    filter(zstart2 == min(zstart2)) %>%
    ungroup()
write_tsv(d, 'lastz-align-10k-gapped-reciprocal-filtered.bed')
```

onto a notebook (`reciprocal_alignment.ipynb`) to check how many bases
are _actually_ unique and how many are gametologs

other things to figure out - how much overlap is
there between the reciprocal alignments?
the reciprocal has less - are these any
alignments in that that aren't in the main alignment?
(although LASTZ should by definition circumvent
the need for reciprocal alignments...)

also check whether the gff coordinates of mt specific genes
(ie fus, mid) overlap with lastz alignments - this
shouldn't be the case! 

also - dotplot for alignment? would be kind of cool...

## 21/12/2018

things to do
- overlap between reciprocal alignments
- GFF coordinates of fus/mid/etc - do these show up in the alignment? (they shouldn't!)
- make a dotplot for the plus-target alignment
    - does this line up with the synteny diagram in de Hoff 2013?
- figure out what's going wrong in the unique base count notebook
    - could create filtered versions of lastz files that concatenate overlapping intervals

### checking mt specific genes

mt+ specific genes include FUS1 and MTA1

mt- specific genes include MIDm and MTD1m

we can check for these in `data/references/final.strict.GFF3` by searching with `less`:

```bash
131574 chromosome_6    phytozome8_0    gene    480316  484847  .       +       .       ID=Cre06.g252750;Name=Cre06.g2
131575 chromosome_6    phytozome8_0    mRNA    480316  484847  .       +       .       ID=PAC:26894229;Name=Cre06.g25
```

let's correct the lastz alignment coordinates in R:

```R
library(tidyverse); library(magrittr)
d <- read_tsv('lastz-align-10k-gapped-filtered.bed')
d %<>%
    mutate(zstart1 = zstart1 + 298299, # GFFs are 1-based
           end1 = end1 + 298299,
           zstart2 = zstart2 + 1,
           end2 = end2 + 1)
```

seems we're clear:

```
> d %>% filter(zstart1 > 470000, end1 < 500000) %>% arrange(zstart1)
# A tibble: 15 x 13
    score        name1 strand1 zstart1   end1   name2 strand2 zstart2   end2
    <int>        <chr>   <chr>   <dbl>  <dbl>   <chr>   <chr>   <dbl>  <dbl>
 1  75172 chromosome_6       +  474234 475088 mtMinus       -  124198 125052
 2  86664 chromosome_6       +  475343 476591 mtMinus       +  305633 306932
 3  13075 chromosome_6       +  477141 477338 mtMinus       -  201014 201219
 4  35043 chromosome_6       +  485701 486169 mtMinus       -  212822 213275
 5  28192 chromosome_6       +  485718 486169 mtMinus       -  161694 162117
 6  50548 chromosome_6       +  489069 489778 mtMinus       +  313533 314289
 7  25975 chromosome_6       +  489136 489775 mtMinus       -  229401 230192
 8  33732 chromosome_6       +  489179 489792 mtMinus       +  150989 151536
 9  11903 chromosome_6       +  489415 489605 mtMinus       +  104620 104812
10  17813 chromosome_6       +  491434 491659 mtMinus       -  171509 171733
11 244590 chromosome_6       +  492109 495281 mtMinus       -  174358 177487
12  34905 chromosome_6       +  495254 495714 mtMinus       -  178222 178675
13  29553 chromosome_6       +  496259 496704 mtMinus       -  180325 180723
14  74262 chromosome_6       +  496757 497657 mtMinus       -  182078 182988
15  32508 chromosome_6       +  497652 498447 mtMinus       -  183795 184325
```

mta1, mid, mtd:

```bash
# MTA1
131808 chromosome_6    phytozome8_0    gene    559654  561075  .       -       .       ID=Cre06.g253000;Name=Cre06.g2
131809 chromosome_6    phytozome8_0    mRNA    559654  561075  .       -       .       ID=PAC:26893117;Name=Cre06.g25

# MID
423581 mtMinus feature gene    286502  287123  .       +       .       gene=MIDm;ness_ID=ADF43190.1;ID=ADF43190.1;Nam
423582 mtMinus feature CDS     286502  286653  .       +       0       ness_ID=ADF43190.1;ID=ADF43190.1.CDS.1;Parent=

# MTD
423600 mtMinus feature gene    308874  311600  .       -       .       gene=MTD1m;ness_ID=ADF43193.1;ID=ADF43193.1;Na
423601 mtMinus feature CDS     308874  309104  .       -       0       ness_ID=ADF43193.1;ID=ADF43193.1.CDS.1;Parent=
```

worth noting in the case of MTA1 that even though the gene is on the rev orientation,
the GFF coordinates are forward based - and so the same coordinates hold in the lastz
file, since it locks the target to the fwd orientation and reports the rev comp query
if necessary

in R:

```R
# mta1
> d %>% filter(zstart1 > 559000, zstart1 < 565000)
# A tibble: 1 x 13
   score        name1 strand1 zstart1   end1   name2 strand2 zstart2   end2
   <int>        <chr>   <chr>   <dbl>  <dbl>   <chr>   <chr>   <dbl>  <dbl>
1 113425 chromosome_6       +  561263 562542 mtMinus       -  287257 288562

# mid
> d %>% filter(zstart2 > 285000, zstart2 < 290000)
# A tibble: 3 x 13
   score        name1 strand1 zstart1   end1   name2 strand2 zstart2   end2
   <int>        <chr>   <chr>   <dbl>  <dbl>   <chr>   <chr>   <dbl>  <dbl>
1  50598 chromosome_6       +  447609 448174 mtMinus       +  288002 288574
2 445456 chromosome_6       +  449521 454907 mtMinus       +  288625 293676
3 113425 chromosome_6       +  561263 562542 mtMinus       -  287257 288562

# mtd
> 345555 - c(308874, 311600) # convert coordinates
[1] 36681 33955

> d %>% filter(zstart2 > 30000, zstart2 < 40000, strand2 == '-')
# A tibble: 4 x 13
  score        name1 strand1 zstart1   end1   name2 strand2 zstart2  end2
  <int>        <chr>   <chr>   <dbl>  <dbl>   <chr>   <chr>   <dbl> <dbl>
1 51697 chromosome_6       +  458413 459135 mtMinus       -   31250 32024
2 14199 chromosome_6       +  592791 592974 mtMinus       -   37480 37663
3 12342 chromosome_6       +  592792 592953 mtMinus       -   39811 39972
4 15822 chromosome_6       +  662469 662690 mtMinus       -   37460 37685
```

verdict: none of the mt-specific genes show up in the lastz alignment! (phew)

### alignment dotplot

from the lastz documentation - the R dotplot is a homemade format
that is engineered to return a dotplot when fed into R's base::plot()

```bash
time ./bin/lastz data/references/mt_plus.fasta data/references/mt_minus.fasta \
--output=data/alignment-lastz/lastz-align-10k-gapped.dotplot \
--hspthresh=10000 \
--format=rdotplot
```

this looks messy - how about higher hsp thresholds?

```bash
time ./bin/lastz data/references/mt_plus.fasta data/references/mt_minus.fasta \
--output=data/alignment-lastz/lastz-align-50k-gapped.dotplot \
--hspthresh=50000 \
--format=rdotplot

time ./bin/lastz data/references/mt_plus.fasta data/references/mt_minus.fasta \
--output=data/alignment-lastz/lastz-align-100k-gapped.dotplot \
--hspthresh=100000 \
--format=rdotplot
```

okay - these plots make clear just how noisy the 10k alignment is - and how much
it doesn't actually line up with the geneious alignment! 

I think we should work with a 50k alignment from here on out. I should have checked
this far sooner - but at least there's an existing codebase to quickly redo these
analyses with. 

```bash
time ./bin/lastz data/references/mt_plus.fasta data/references/mt_minus.fasta \
--output=data/alignment-lastz/lastz-align-50k-gapped.bed \
--hspthresh=50000 \
--format=general

time ./bin/lastz data/references/mt_plus.fasta data/references/mt_minus.fasta \
--output=data/alignment-lastz/lastz-align-50k-gapped.maf \
--hspthresh=50000 \
--format=maf

time bash main.sh
```

## 25/12/2018

did I use the wrong VCFs again?...

I've been using the VCFs marked `gGVCFs` and not `HC`, but
it appears that `HC` is the right one? 

let's try running these things with the other VCFs in the meantime
at least:


```bash
cd data/references
ln -sv ../../j_old/analysis/map_MT-specific/MT-specific_VCFs/all_quebec.mtMinus.HC.vcf.gz* .
ln -sv ../../j_old/analysis/map_MT-specific/MT-specific_VCFs/all_quebec.mtPlus.HC.vcf.gz* .

time ./bin/vcf2fasta.py -v data/references/all_quebec.mtPlus.HC.vcf.gz \
-r data/references/mtPlus_ref.chromosome_6.fasta \
-i chromosome_6:298299-826737 \
-s CC2936 CC2937 CC3060 CC3064 CC3065 CC3068 CC3071 CC3076 CC3086 \
--min_GQ 30 > data/aligned-fastas/plus_strains_ref_HC.fasta

time ./bin/vcf2fasta.py -v data/references/all_quebec.mtMinus.HC.vcf.gz \
-r data/references/mtMinus_ref.chromosome_6_and_mtMinus.fasta \
-i mtMinus:1-345555 \
-s CC2935 CC2938 CC3059 CC3061 CC3062 CC3063 CC3073 CC3075 CC3079 CC3084 \
--min_GQ 30 > data/aligned-fastas/minus_strains_ref_HC.fasta
```

holy moly - look at all this data we were missing in the plus!

```python
>>> from Bio import SeqIO
>>> fname = 'data/aligned-fastas/plus_strains_ref_HC.fasta'
>>> d = [s for s in SeqIO.parse(fname, 'fasta')]
>>> [len(s) for s in d]
[528439, 528439, 528439, 528439, 528439, 528439, 528439, 528439, 528439]
>>> [str(s.seq).count('N') for s in d]
[28498, 28379, 28852, 29052, 28842, 28493, 28841, 28613, 28038]
>>> old_fname = 'data/aligned-fastas/plus_strains_ref.fasta'
>>> d2 = [s for s in SeqIO.parse(old_fname, 'fasta')]
>>> [str(s.seq).count('N') for s in d2]
[191768, 187475, 197104, 205915, 203036, 184006, 202624, 195931, 167198]
```

and in the minus:

```python
>>> from Bio import SeqIO
>>> fname = 'data/aligned-fastas/minus_strains_ref_HC.fasta'
>>> old_fname = 'data/aligned-fastas/minus_strains_ref.fasta'
>>> s = [s for s in SeqIO.parse(fname, 'fasta')]
>>> d = [s for s in SeqIO.parse(fname, 'fasta')]
>>> d2 = [s for s in SeqIO.parse(old_fname, 'fasta')]
>>> print([len(s) for s in d], [len(s) for s in d2])
[345555, 345555, 345555, 345555, 345555, 345555, 345555, 345555, 345555, 345555] [345555, 345555, 345555, 345555, 345555, 345555, 345555, 345555, 345555, 345555]
>>> [str(s.seq).count('N') for s in d]
[3346, 3225, 3440, 3631, 3292, 3711, 3306, 3426, 3480, 3510]
>>> [str(s.seq).count('N') for s in d2]
[36370, 35274, 42615, 41411, 42410, 47596, 40695, 43703, 44428, 46323]
```

we'll create a modified version of main called `main2.sh` for
now, incorporating the new filenames and making sure not
to overwrite old ones - let's see if this makes a difference

## 26/12/2018

it didn't - the rho values are still being skewed.

we should get the DP values across the mt locus and see
what these look like. here's a quick script to do just that:

```python
import vcf
import sys
from tqdm import tqdm

fname = sys.argv[-2]
outname = sys.argv[-1]

with open(outname, 'w') as f:
    first_iteration = True
    records = vcf.Reader(filename = fname, compressed = True)
    region = records.fetch('chromosome_6', 298000, 830000) # mt plus
    for record in tqdm(region):
        if first_iteration:
            sample_names = [c.sample for c in record.samples]
            f.write('POS ' + ' '.join(sample_names) + '\n')
            col_count = len(sample_names)
            first_iteration = False
        try:
            depth_vals = [sample['DP'] for sample in record.samples]
        except AttributeError:
            continue
        for i, val in enumerate(depth_vals):
            if isinstance(val, int):
                depth_vals[i] = str(depth_vals[i])
            elif val is None:
                depth_vals[i] = 'NA'
        assert len(depth_vals) == col_count
        position = str(record.POS)
        f.write(position + ' '  + ' '.join(depth_vals) + '\n')
```            

```bash
time python3.5 depth_check.py \
data/references/all_quebec.mtPlus.HC.vcf.gz \
depth_test.txt

time python3.5 depth_check.py \
data/references/all_quebec.mtPlus.gGVCFs.vcf.gz \
depth_test_gvcf.txt
```

let's try redoing everything _again_, but this time applying
a min DP of 30. this is a reasonably high filter but
let's see if that does away with the weirdness we're seeing.
onto `main3.sh`...

update - same culprits - nothing's changed.

wait a second - de Hoff mentioned autosomal translocations
in the mt+ hap, and the depth for that region seems
much higher than surrounding (plotted in RStudio locally) -
what do we get from a lastz dotplot between that chromosome
and the mt+?

the de Hoff paper says it was chromosome 16 - and so:

```bash
time ./bin/lastz data/references/mt_plus.fasta \
../genomewide_recombination/data/fastas/reference/chromosome_16.fasta \
--output=autosomal_test.dotplot \
--hspthresh=10000 \
--format=rdotplot

time ./bin/lastz data/references/mt_plus.fasta \
../genomewide_recombination/data/fastas/reference/chromosome_16.fasta \
--output=autosomal_test.bed \
--hspthresh=10000 \
--format=general

# maf - high thresh
time ./bin/lastz data/references/mt_plus.fasta \
../genomewide_recombination/data/fastas/reference/chromosome_16.fasta \
--output=autosomal_test.maf \
--hspthresh=100000 \
--format=maf

```

see also fig. 2 in de hoff

overall - it might seem autosomal translocations might help
explain the huge variance in sequencing depth, thus potentially
affecting our results - but the two primary regions
that seemed to have high recombination (535k and 424k) seem to be
just upstream of the actual translocation regions (something
my own lastz alignment above confirms). what could be
going on here?...

## 27/12/2018

another thing to consider - variance in SNP density

if SNP density is higher in these regions, we'll likely see higher recombination
rate estimates

we can calculate SNP density in windows with a quick python script:

```python
import vcf
import sys
from tqdm import tqdm

fname = sys.argv[-3]
outname = sys.argv[-2]
windowsize = sys.argv[-1]

with open(outname, 'w') as f:
    first_iteration = True
    for segment in tqdm(range(298000, 830000, windowsize)):
        if first_iteration:
            f.write('start end snp_count\n')
            first_iteration = False

        start = segment
        end = segment + windowsize

        v = vcf.Reader(filename = fname, compressed = True)
        region = v.fetch('chromosome_6', start, end)

        snps = len([record.is_snp for record in region])
        f.write(' '.join([start, end, snps] + '\n')

```

followed by

```bash
cd autosomal_test
time python3.5 snp_density.py \
../data/references/all_quebec.mtPlus.HC.vcf.gz \
snp_density_out.txt \
2000
```

this doesn't seem to be the biggest predictor:

```R
> d %>% arrange(desc(snp_count)) %>% head(20)
# A tibble: 20 x 3
    start    end snp_count
    <int>  <int>     <int>
 1 662000 664000       276
 2 320000 322000       275
 3 328000 330000       239
 4 338000 340000       233
 5 332000 334000       218
 6 540000 542000       213
 7 814000 816000       212
 8 560000 562000       209
 9 306000 308000       199
10 322000 324000       199
11 536000 538000       198
12 312000 314000       194
13 360000 362000       192
14 596000 598000       191
15 348000 350000       189
16 358000 360000       187
17 334000 336000       185
18 354000 356000       184
19 346000 348000       183
20 574000 576000       182
```

at a 10kb windowsize, that region definitely seems to have really high SNP density -
but so does the region from 298k to 400k (ie the T domain iirc - which I would
expect to have somewhat higher SNP density)

one might hope that the region in the T domain has higher density due to recombination
while the nongametologous region is larger due to the translocation - but is there
any way to really be sure? 


## 28/12/2018

why don't we see the C domain in the lastz alignment?

de Hoff claims that the C domain encompasses 116 kb of sequence at the end of the mt+ locus,
and is syntenic with the mt- locus, but we don't see alignments there.

what if we aligned the final 150 kb of each locus vs one another?

```python
from Bio import SeqIO
plus = SeqIO.read('data/references/mt_plus.fasta', 'fasta')
minus = SeqIO.read('data/references/mt_minus.fasta', 'fasta')
plus_final_150 = plus[len(plus) - 150000:]
minus_final_150 = minus[len(minus) - 150000:]
with open('plus_final_150.fasta', 'w') as f:
    SeqIO.write(plus_final_150, f, 'fasta')
with open('minus_final_150.fasta, 'w') as f:
    SeqIO.write(minus_final_150, f, 'fasta')
```

```bash
mkdir c-domain
mv -v *150.fasta c-domain/
./bin/lastz c-domain/plus_final_150.fasta c-domain/minus_final_150.fasta \
--output=c-domain/alignment_10k.bed \
--hspthresh=10000 \
--format=general
```

after playing around, it seems a 30k score does it right.

the key alignment is at the end here:


```bash
 15 1389659 chromosome_6    +       150000  104633  119541  mtMinus +       150000  130541  145437  14658/14827
 16 204537  chromosome_6    +       150000  140345  142902  mtMinus +       150000  119072  121435  2264/2346
 17 165936  chromosome_6    +       150000  141121  142911  mtMinus +       150000  0       1789    1751/1784
 18 399343  chromosome_6    +       150000  145338  150000  mtMinus +       150000  145381  150000  4404/4553
```

do we... do have to do the entire mt locus alignment again...

```bash
time ./bin/lastz data/references/mt_plus.fasta data/references/mt_minus.fasta \
--output=data/alignment-lastz/lastz-align-30k-gapped.dotplot \
--hspthresh=30000 \
--format=rdotplot

time ./bin/lastz data/references/mt_plus.fasta data/references/mt_minus.fasta \
--output=data/alignment-lastz/lastz-align-30k-gapped.bed \
--hspthresh=30000 \
--format=general

time bash main.sh
```

## 31/12/2018

given the DP investigations I did locally, let's try filtering
sites that are over 3x the mean DP to get rid of any
potential paralogous regions that have mapped to the mt.

vcf2fasta doesn't have a flag for this, but we
can filter using a quick python script. this script
should

1. calculate mean DP across the mt locus
2. iterate through the VCF again and write records where DP < 3 * mean

after that, we can use vcf2fasta to filter, and get
a variant of `main.sh` running again.

```python
import vcf
import sys

allele = sys.argv[-4]
fname = sys.argv[-3]
outname = sys.argv[-2]
depth_multiplier = sys.argv[-1]

if allele == 'plus':
    chrom = 'chromosome_6'
    start = 298299
    end = 826737
elif allele == 'minus':
    chrom = 'mtMinus'
    start = 1
    end = 345555

# get mean DP
dp_vals = []
vcfin = vcf.Reader(filename = fname, compressed = True)
dp_vals = [record.INFO['DP'] for record
           in vcfin.fetch(chrom, start, end)]
mean_dp = sum(dp_vals) / len(dp_vals)

dp_threshold = mean_dp * depth_multiplier

# write new VCF
written = 0
counter = 0
with open(outfile, 'w') as f:
    vcfin = vcf.Reader(filename = fname, compressed = True)
    writer = vcf.Writer(f, vcfin)
    for record in vcfin.fetch(chrom, start, end):
        counter += 1
        if record.INFO['DP'] < dp_threshold:
            written += 1
            writer.write_record(record)
    
print('Completed filtering.')
print('Records below', dp_threshold, 'were removed.')
print(counter, 'records iterated over.')
print(written, 'records passed filtering.')
```
    
```bash
python3.5 analysis/alignment-lastz/filter_max_dp.py \
plus data/references/all_quebec.mtPlus.HC.vcf.gz \
data/references/all_quebec.mtPlus.HC.DP3x.vcf.gz 3
```

so the mean DP is driven too far upward already - this
needs to be done on a strain-by-strain basis.

```bash
# new and improved!
time python3.5 analysis/alignment-lastz/filter_max_dp.py \
--filename data/references/all_quebec.mtPlus.HC.vcf.gz \
--outfile data/references/all_quebec.mtPlus.HC.DP2x.vcf \
--allele plus \
--dp_multiplier 2

time python3.5 analysis/alignment-lastz/filter_max_dp.py \
--filename data/references/all_quebec.mtMinus.HC.vcf.gz \
--outfile data/references/all_quebec.mtMinus.HC.DP2x.vcf \
--allele minus \
--dp_multiplier 2
```

```bash
bgzip data/references/all_quebec.mtPlus.HC.DP2x.vcf
bgzip data/references/all_quebec.mtMinus.HC.DP2x.vcf
tabix -p vcf data/references/all_quebec.mtPlus.HC.DP2x.vcf.gz
tabix -p vcf data/references/all_quebec.mtMinus.HC.DP2x.vcf.gz
```

and now we'll have another stab at a version of `main.sh`,
which also creates new fastas:

```bash
time bash main_dp.sh
```
