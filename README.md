# tRNA-m5C
Scripts for tRNA BS-seq alignment.

Requirement:
pysam (0.15.2 tested)

biopython

numpy

scipy

### Metadata

(1) mature tRNA sequences (from GtRNAdb)

(2) mRNA reference sequences (from Ensembl/UCSC/etc.)

(3) rRNA reference sequences (from SILVA)

### Reference preparation

Note: this step aims to remove duplicate sequences from the mature tRNA sequence file from GtRNAdb. For example, tRNA-Arg-ACG-1-1 and tRNA-Arg-ACG-1-2 are copies of tRNA-Arg-ACG-1 and have an identical sequence, so they should be considered as ONE sequence in the analysis. Then the CCA tail is appended to the mature tRNA sequences to ensure accurate alignment.

python format_mature_fasta.py <mature.tRNA.fa> > <mature.tRNA.format.fa>

\# Do it yourself, seperate mitochondrial tRNA sequence and mRNA sequence from Ensembl cDNA fasta file, add CCA tail for them.

\# Suppose you have <MT.tRNA.fa> and <mRNA.fa>

cat <mature.tRNA.format.fa> <MT.tRNA.fa> > <tRNA.fa>

python add_CCA.py <tRNA.fa> > <tRNA.CCA.fa>

cat <tRNA.CCA.fa> <mRNA.fa> <rRNA.fa> > <reference.fa>

python fasta_c2t.py -i <reference.fa> > <reference.c2t.fa>

\# Build index

bowtie2-build <reference.c2t.fa>

### The pipeline

#### 1. Trim adapter, we use tRNA kit from Vazyme: universial_adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", small_RNA_adapter = "GATCGTCGGACTGTAGAACTCTGAAC"

cutadapt --max-n 1 -m 18 -e 0.25 --trim-n -q 20 --trimmed-only -a {universial_adapter} -A {small_RNA_adapter} -o <read1.cutadapt.fastq> -p <read2.cutadapt.fastq> <read1.fastq> <read2.fastq>

#### 2. C2T/G2A conversion

fastq_c2a.py -i <read1.cutadapt.fastq> > <read1.cutadapt.c2t.fastq>

fastq_g2a.py -i <read2.cutadapt.fastq> > <read2.cutadapt.g2a.fastq>

#### 3. Alignment. -X 80 for maximum PE alignment distance; -k 50 for a shorter running time, you can change it to -a.
bowtie2 -x {bowtie2_index} --end-to-end --no-mixed --norc --no-unal -k 50 -p 20 -S <bowtie2.sam> -1 <read1.cutadapt.c2t.fastq> -2 <read2.cutadapt.g2a.fastq>

#### 4. Rescue Cs in the SAM

python tRNA_bam_recovery_PE.py -i <bowtie2.sam> -o <bowtie2.bam> -f <read1.cutadapt.fastq> -r <read2.cutadapt.fastq>

#### 5. Get rid of multiple alignments

python filter_tRNA_alignments.v2.py -i <tRNA.bam> -o <tRNA.filtered.bam>

#### 6. Sort and index, adjust -@ and -m for a better performence

samtools sort -@ 4 -m 4G -o <tRNA.filtered.sorted.bam><tRNA.filtered.bam>

samtools index <tRNA.filtered.sorted.bam>

#### 7. Pileup and enjoy the result

python pileup_tRNA_bases_PE.py -b <tRNA.filtered.sorted.bam> -r <tRNA.CCA.fa> -o <pileup.res.txt>










