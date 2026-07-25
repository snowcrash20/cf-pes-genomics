#!/usr/bin/env python3

import csv
from pathlib import Path


RESULT_DIR = Path("results/representative_real")
FORBIDDEN_TEXT = (
    "161_AL",
    "214_AL",
    "MAC_Aer",
    "k141_",
    "Homo sapiens",
    "9606",
)


def read_csv(name: str) -> list[dict[str, str]]:
    with (RESULT_DIR / name).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


compact = read_csv("compact_contig_taxonomy.csv")
reports = read_csv("kraken_report_summary.csv")
species = read_csv("top_bacterial_species.csv")
features = read_csv("P749_reference_features_of_interest.csv")

all_text = "\n".join(
    path.read_text(encoding="utf-8", errors="strict")
    for path in RESULT_DIR.iterdir()
    if path.suffix in {".csv", ".md"}
)

for forbidden in FORBIDDEN_TEXT:
    if forbidden in all_text:
        raise AssertionError(f"Forbidden source identifier found: {forbidden}")

totals = {}
for sample in ("Onset", "Resolution"):
    totals[sample] = sum(
        int(row["contig_count"]) for row in compact if row["sample"] == sample
    )

assert totals == {"Onset": 812, "Resolution": 2481}
assert sum(
    int(row["contig_count"])
    for row in compact
    if row["taxon_name"] == "Host-classified"
) == 36
assert "contig_id" not in compact[0]
assert "hit_list" not in compact[0]
assert len(species) == 10
assert any(row["gene"] == "gyrA" for row in features)
assert any(row["gene"] == "mexA" for row in features)
assert any(row["taxon_name"] == "Host-classified" for row in reports)

figure = RESULT_DIR / "representative_kraken_species.png"
assert figure.exists() and figure.stat().st_size > 10_000

print("Representative outputs are compact and deidentified.")
