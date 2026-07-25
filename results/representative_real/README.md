# Representative real results

> **Approval required before public release:** these compact tables were derived
> from two patient-associated project samples. Confirm public sharing with the
> project supervisor or data owner before pushing this directory to GitHub.

This directory demonstrates the real structure and biological range of the
project outputs without publishing raw reads, assembled sequences, contig
identifiers, per-contig k-mer hit lists, or original study sample IDs.

Files:

- `compact_contig_taxonomy.csv`: per-sample aggregate of exact Kraken2 contig
  assignments, with counts and cumulative sequence length.
- `kraken_report_summary.csv`: compact long-format representation of the
  Kraken2 reports.
- `top_bacterial_species.csv`: the five most abundant bacterial species across
  the two representative samples.
- `representative_kraken_species.png`: portfolio figure comparing the
  deidentified onset and resolution samples.
- `P749_reference_features_of_interest.csv`: small subset of the later Bakta
  annotation containing recurrent features discussed in the project.

Host-classified assignments are collapsed to `Host-classified`. No host
sequences are included. The original `.out` files, full Bakta TSV, FASTA files,
and original sample identifiers are deliberately excluded.

The Bakta subset came from a later annotation generated with Bakta 1.11.3 and
database 6.0. It is therefore labelled as a later reference annotation and
should not be presented as the exact Bakta 1.7.0 annotation used during the
original analysis.

These files can be regenerated from private source files with:

```bash
python scripts/prepare_representative_results.py \
  --onset-out /private/path/onset.out \
  --onset-report /private/path/onset.report \
  --resolution-out /private/path/resolution.out \
  --resolution-report /private/path/resolution.report \
  --bakta-tsv /private/path/P749_bakta.tsv \
  --out-dir results/representative_real
```
