# PCR Duplicates

### Why Remove Them?

Differential expression

Genome assembly

### How to Identify them?

**SAM Format** - same position, same chromosome, same strand

**Soft Clipping** - Definite not matching at the ends of a read

Cigar String 

- M - Match or mismatch
- S - Soft Clipping

Match position starts at M, not S for something like 2S12M

[zenfractal.com](http://zenfractal.com/2013/06/19/playing-with-matches/) "Playing with Matches"

**UMI** - Unique Molecular Index, PCR duplicates have the same UMI

​	We could have actual duplicate molecules that aren't PCR duplicates, UMIs differentiate them

#### Tools

* Picard

* Umitools

For Paired-end data, if a pairs of reads have the same UMIs and same position it's likely a PCR duplicate

Never assume UMIs are all at even levels

If you've quality-trimmed the UMI (5') end that will affect start position 

## Our Algorithm

Given a SAM file of uniquely-mapped reads, remove all PCR duplicates

Use Samtools sort

Adjust for soft clipping

Start with single-end, expand to paired-end

Start with known UMIs

Millions of reads, don't just load everything into memory

### Part 1

Define the problem

Write Psuedocode

Test examples

Return statement



