# Data layout

Only synthetic demonstration inputs are version-controlled.

For real analysis, use the following local layout:

```text
data/
├── reference/
│   ├── P749.fasta
│   ├── P749_bakta.gbk
│   └── P749_bakta.gff3
├── snp_tables/
│   └── <sample>_snps_breseq.tsv
├── megahit/
│   └── <sample>/final.contigs.fa
└── kraken2/
    ├── <sample>.out
    └── <sample>.report
```

The `.gitignore` excludes these real-data directories and common sequencing
formats. Do not upload patient-derived reads, intermediate alignments, or
controlled reference assets unless their redistribution is explicitly
permitted.
