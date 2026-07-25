#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/.." && pwd)"

cd "${repo_dir}"

Rscript R/snp_frequency_analysis.R \
  --snp-dir data/example/snp_tables \
  --gff data/example/P749_subset.gff3 \
  --metadata config/example_samples.tsv \
  --aliases config/product_aliases.tsv \
  --out-dir results/example/snp \
  --qual-threshold 10 \
  --fixed-threshold 0.95 \
  --top-n 3

Rscript R/unmapped_contig_taxonomy.R \
  --megahit-dir data/example/megahit \
  --kraken-dir data/example/kraken2 \
  --metadata config/example_samples.tsv \
  --out-dir results/example/contigs \
  --top-n 3

test -s results/example/snp/tables/annotated_snps.csv
test -s results/example/snp/tables/filter_summary.csv
test -s results/example/snp/figures/top_features_per_sample.png
test -s results/example/contigs/tables/contig_taxonomy.csv
test -s results/example/contigs/figures/unmapped_contig_taxonomy.png

Rscript tests/validate_example_outputs.R

echo "Example workflows completed successfully."
