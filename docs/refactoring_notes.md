# Refactoring notes

The uploaded scripts captured the full exploratory analysis but were tied to
one Windows workstation and one fixed set of sample names. The repository
versions preserve the analytical intent while changing the interface.

| Original pattern                                     | Repository implementation                             |
| ---------------------------------------------------- | ----------------------------------------------------- |
| `setwd()` and `~/Rstudio/...` paths                  | Required command-line paths                           |
| Figures printed interactively                        | Named PNG/PDF files written to an output directory    |
| Patient pairs written one by one                     | Pairwise comparisons generated from metadata          |
| Only overlapping SNPs retained                       | All SNPs retained; non-overlaps labelled `Intergenic` |
| Both `gene` and `CDS` overlaps could duplicate a SNP | CDS preferred deterministically                       |
| Undefined intermediate `snps_all_labeled`            | Summaries derived directly from validated data        |
| Fixed sample-name lookup inside code                 | Reusable tab-delimited metadata file                  |
| Kraken rows joined to FASTA lengths                  | FASTA contigs used as the left side of the join       |
| Missing Kraken classifications could disappear       | Missing classifications retained explicitly           |
| No executable public example                         | Synthetic inputs plus a CI smoke test                 |

The original plotting-specific colour choices and report-only `gt` tables were
not retained in the core workflows. The public scripts prioritise stable CSV
outputs and portable figures that can be restyled downstream.
