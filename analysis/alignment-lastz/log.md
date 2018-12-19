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


