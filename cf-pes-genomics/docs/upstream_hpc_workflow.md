# Upstream HPC workflow

This document records the command-line stages that produce the stable inputs
consumed by the R scripts. Adapt scheduler directives, module names, paths, and
CPU or memory requests to the local cluster.

## 1. Read quality control and trimming

Initial and post-trimming quality were assessed with FastQC. The original
analysis used Trimmomatic in paired-end mode with:

```bash
trimmomatic PE \
  -threads "${SLURM_CPUS_PER_TASK}" \
  sample_R1.fastq.gz sample_R2.fastq.gz \
  sample_R1_paired.fastq.gz sample_R1_unpaired.fastq.gz \
  sample_R2_paired.fastq.gz sample_R2_unpaired.fastq.gz \
  ILLUMINACLIP:TruSeq3-PE-2.fa:2:20:10:2:true \
  LEADING:5 \
  TRAILING:5 \
  SLIDINGWINDOW:4:20 \
  MINLEN:50
```

Keep FastQC HTML and ZIP outputs as provenance, but do not commit large raw-read
files to Git.

## 2. Polymorphism-aware variant calling

breseq was run in polymorphism mode against an annotated P749 reference:

```bash
breseq \
  -p \
  -j "${SLURM_CPUS_PER_TASK}" \
  -r P749_bakta.gbk \
  -o sample_breseq \
  sample_R1_paired.fastq.gz \
  sample_R2_paired.fastq.gz
```

Use `breseq --help` for the installed release if a long option differs. The
project used breseq 0.36.1.

Export SNP records from the resulting VCF to one tab-delimited file per sample.
The R interface deliberately starts at the following eight-column contract:

```text
CHROM  POS  REF  ALT  QUAL  AF  AD  DP
```

bcftools query syntax depends on whether a VCF stores `AF`, `AD`, and `DP` in
INFO or FORMAT fields. Inspect the VCF header first:

```bash
bcftools view --header-only sample.vcf | less
```

For per-sample FORMAT fields, the export resembles:

```bash
bcftools view --types snps sample.vcf |
  bcftools query \
    --format '%CHROM\t%POS\t%REF\t%ALT\t%QUAL[\t%AF\t%AD\t%DP]\n' \
    > sample_snps_breseq.tsv
```

Confirm that the output contains exactly eight columns before running the R
analysis. If the tags are INFO fields, use `%INFO/AF`, `%INFO/AD`, and
`%INFO/DP` instead.

## 3. Retrieve unmapped paired reads

Build the P749 Bowtie2 index once, then retain pairs that fail to align:

```bash
bowtie2-build P749.fasta P749

bowtie2 \
  --very-sensitive \
  -p "${SLURM_CPUS_PER_TASK}" \
  -x P749 \
  -1 sample_R1_paired.fastq.gz \
  -2 sample_R2_paired.fastq.gz \
  --un-conc-gz sample_unmapped.fastq.gz \
  --sam-no-qname-trunc \
  --sam-append-comment \
  -S sample.sam
```

Bowtie2 writes the paired unmapped files using its numbered filename expansion.
Record any sample-specific relaxation of mismatch or score thresholds. Such a
change means that the resulting unmapped-read sets are not directly equivalent
unless the exception is documented.

## 4. Assemble unmapped reads

The original analysis used MEGAHIT 1.2.9, an explicit k-mer list, and a minimum
contig length of 1,000 bp:

```bash
megahit \
  -t "${SLURM_CPUS_PER_TASK}" \
  -m 0.9 \
  --k-list 21,29,39,59,79,99,119 \
  --min-contig-len 1000 \
  -1 sample_unmapped.fastq.1.gz \
  -2 sample_unmapped.fastq.2.gz \
  -o megahit/sample
```

Use either an explicit k-mer list or a preset after checking the installed
MEGAHIT help. Presets can override individual k-mer settings in some releases.

## 5. Classify assembled contigs

Run Kraken2 against a documented database and retain both per-contig output and
the sample report:

```bash
kraken2 \
  --db "${KRAKEN2_DB}" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --use-names \
  --report kraken2/sample.report \
  --output kraken2/sample.out \
  megahit/sample/final.contigs.fa
```

For reproducibility, save at least:

- Kraken2 version;
- database type and build date;
- included reference libraries;
- database checksum or release identifier; and
- command-line options, especially any confidence threshold.

The database itself is large and should normally be referenced by provenance
rather than committed to this repository.
