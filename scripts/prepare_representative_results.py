#!/usr/bin/env python3

"""Create compact, deidentified portfolio outputs from private Kraken/Bakta files."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


csv.field_size_limit(10_000_000)

HOST_TAXIDS = {"9605", "9606"}
FEATURE_GENES = {
    "alg8",
    "algE",
    "fusA",
    "gyrA",
    "hcpA",
    "mexA",
    "morA",
    "ssuD",
}
FEATURE_PRODUCTS = (
    "DNA-Binding protein G5P",
    "Ig-like domain-containing protein",
    "LysR family transcriptional regulator",
    "hemagglutination protein",
    "hemagglutinin",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate private per-contig Kraken2 outputs into compact, "
            "deidentified portfolio tables and a comparison figure."
        )
    )
    parser.add_argument("--onset-out", required=True, type=Path)
    parser.add_argument("--onset-report", required=True, type=Path)
    parser.add_argument("--resolution-out", required=True, type=Path)
    parser.add_argument("--resolution-report", required=True, type=Path)
    parser.add_argument("--bakta-tsv", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--top-n", default=5, type=int)
    return parser.parse_args()


def sanitize_taxon(taxid: str, taxon_name: str) -> tuple[str, str]:
    if taxid in HOST_TAXIDS or taxon_name.lower().startswith("homo"):
        return "", "Host-classified"
    return taxid, taxon_name


def parse_taxon_field(value: str) -> tuple[str, str]:
    match = re.match(r"^(.*?)\s+\(taxid\s+([0-9]+)\)$", value)
    if match:
        taxon_name, taxid = match.groups()
        return sanitize_taxon(taxid, taxon_name.strip())
    return sanitize_taxon(value, "")


def aggregate_kraken_out(
    path: Path, sample: str
) -> tuple[list[dict[str, object]], int]:
    aggregate: dict[tuple[str, str, str], dict[str, int]] = defaultdict(
        lambda: {"contig_count": 0, "total_bp": 0}
    )
    total_contigs = 0
    total_bp = 0

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if len(row) != 5:
                raise ValueError(
                    f"{path.name}:{line_number} has {len(row)} columns; expected 5."
                )

            classification, _, taxon_field, length_text, _ = row
            length_bp = int(length_text)
            taxid, taxon_name = parse_taxon_field(taxon_field)
            if classification == "U":
                taxid, taxon_name = "0", "Unclassified"

            key = (classification, taxid, taxon_name)
            aggregate[key]["contig_count"] += 1
            aggregate[key]["total_bp"] += length_bp
            total_contigs += 1
            total_bp += length_bp

    rows: list[dict[str, object]] = []
    for (classification, taxid, taxon_name), values in aggregate.items():
        rows.append(
            {
                "sample": sample,
                "classification": classification,
                "taxid": taxid,
                "taxon_name": taxon_name,
                "contig_count": values["contig_count"],
                "total_bp": values["total_bp"],
                "percent_contigs": round(
                    100 * values["contig_count"] / total_contigs, 4
                ),
                "percent_bp": round(100 * values["total_bp"] / total_bp, 4),
            }
        )

    rows.sort(key=lambda row: (-int(row["total_bp"]), str(row["taxon_name"])))
    return rows, total_contigs


def parse_kraken_report(path: Path, sample: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if len(row) != 6:
                raise ValueError(
                    f"{path.name}:{line_number} has {len(row)} columns; expected 6."
                )

            percent, clade_count, direct_count, rank, taxid, taxon_name = row
            taxid, taxon_name = sanitize_taxon(taxid, taxon_name.strip())
            rows.append(
                {
                    "sample": sample,
                    "percent": float(percent),
                    "clade_contigs": int(clade_count),
                    "direct_contigs": int(direct_count),
                    "rank": rank,
                    "taxid": taxid,
                    "taxon_name": taxon_name,
                }
            )

    return rows


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def select_top_species(
    report_rows: list[dict[str, object]],
    sample_totals: dict[str, int],
    top_n: int,
) -> list[dict[str, object]]:
    species_rows = [
        row
        for row in report_rows
        if row["rank"] == "S"
        and row["taxon_name"] not in {"Host-classified", "Unclassified"}
    ]

    totals: dict[str, int] = defaultdict(int)
    by_sample: dict[tuple[str, str], int] = defaultdict(int)
    for row in species_rows:
        name = str(row["taxon_name"])
        sample = str(row["sample"])
        count = int(row["clade_contigs"])
        totals[name] += count
        by_sample[(sample, name)] += count

    top_species = [
        name
        for name, _ in sorted(totals.items(), key=lambda item: (-item[1], item[0]))[
            :top_n
        ]
    ]

    compact: list[dict[str, object]] = []
    for taxon_name in top_species:
        for sample in ("Onset", "Resolution"):
            count = by_sample[(sample, taxon_name)]
            compact.append(
                {
                    "sample": sample,
                    "taxon_name": taxon_name,
                    "contig_count": count,
                    "percent_of_sample_contigs": round(
                        100 * count / sample_totals[sample], 4
                    ),
                }
            )

    return compact


def plot_top_species(path: Path, rows: list[dict[str, object]]) -> None:
    taxa = list(dict.fromkeys(str(row["taxon_name"]) for row in rows))
    onset = {
        str(row["taxon_name"]): float(row["percent_of_sample_contigs"])
        for row in rows
        if row["sample"] == "Onset"
    }
    resolution = {
        str(row["taxon_name"]): float(row["percent_of_sample_contigs"])
        for row in rows
        if row["sample"] == "Resolution"
    }

    y_positions = list(range(len(taxa)))
    height = 0.36
    figure, axis = plt.subplots(figsize=(10, 5.8))
    axis.barh(
        [position + height / 2 for position in y_positions],
        [onset.get(taxon, 0) for taxon in taxa],
        height=height,
        label="Onset",
        color="#2878B5",
    )
    axis.barh(
        [position - height / 2 for position in y_positions],
        [resolution.get(taxon, 0) for taxon in taxa],
        height=height,
        label="Resolution",
        color="#D95F02",
    )

    axis.set_yticks(y_positions)
    axis.set_yticklabels(taxa, fontstyle="italic")
    axis.invert_yaxis()
    axis.set_xlabel("Percentage of assembled contigs")
    figure.suptitle(
        "Representative bacterial species assignments",
        x=0.01,
        y=0.98,
        horizontalalignment="left",
        fontsize=15,
        fontweight="bold",
    )
    figure.text(
        0.01,
        0.925,
        "Top species across two deidentified longitudinal samples",
        fontsize=10,
        color="#555555",
    )
    axis.grid(axis="x", color="#D9D9D9", linewidth=0.8)
    axis.set_axisbelow(True)
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.legend(frameon=False, loc="lower right")
    figure.text(
        0.01,
        0.01,
        "Host-classified contigs are excluded; counts are derived from Kraken2 reports.",
        fontsize=8.5,
        color="#555555",
    )
    figure.tight_layout(rect=(0, 0.05, 1, 0.89))
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(figure)


def extract_reference_features(path: Path) -> list[dict[str, str]]:
    header: list[str] | None = None
    rows: list[dict[str, str]] = []

    with path.open("r", encoding="utf-8", newline="") as handle:
        for line in handle:
            if line.startswith("#Sequence Id"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if line.startswith("#"):
                continue
            if header is None:
                raise ValueError("Bakta TSV header was not found.")

            values = line.rstrip("\n").split("\t")
            row = dict(zip(header, values, strict=True))
            gene = row["Gene"]
            product = row["Product"]
            if gene in FEATURE_GENES or any(
                phrase.lower() in product.lower() for phrase in FEATURE_PRODUCTS
            ):
                rows.append(
                    {
                        "sequence_id": row["Sequence Id"],
                        "start": row["Start"],
                        "stop": row["Stop"],
                        "strand": row["Strand"],
                        "locus_tag": row["Locus Tag"],
                        "gene": gene,
                        "product": product,
                    }
                )

    return rows


def main() -> None:
    args = parse_args()
    if args.top_n < 1:
        raise ValueError("--top-n must be at least 1.")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    onset_compact, onset_total = aggregate_kraken_out(args.onset_out, "Onset")
    resolution_compact, resolution_total = aggregate_kraken_out(
        args.resolution_out, "Resolution"
    )
    compact_out = onset_compact + resolution_compact
    write_csv(
        args.out_dir / "compact_contig_taxonomy.csv",
        compact_out,
        [
            "sample",
            "classification",
            "taxid",
            "taxon_name",
            "contig_count",
            "total_bp",
            "percent_contigs",
            "percent_bp",
        ],
    )

    report_rows = parse_kraken_report(args.onset_report, "Onset")
    report_rows += parse_kraken_report(args.resolution_report, "Resolution")
    write_csv(
        args.out_dir / "kraken_report_summary.csv",
        report_rows,
        [
            "sample",
            "percent",
            "clade_contigs",
            "direct_contigs",
            "rank",
            "taxid",
            "taxon_name",
        ],
    )

    top_species = select_top_species(
        report_rows,
        {"Onset": onset_total, "Resolution": resolution_total},
        args.top_n,
    )
    write_csv(
        args.out_dir / "top_bacterial_species.csv",
        top_species,
        [
            "sample",
            "taxon_name",
            "contig_count",
            "percent_of_sample_contigs",
        ],
    )
    plot_top_species(args.out_dir / "representative_kraken_species.png", top_species)

    reference_features = extract_reference_features(args.bakta_tsv)
    write_csv(
        args.out_dir / "P749_reference_features_of_interest.csv",
        reference_features,
        [
            "sequence_id",
            "start",
            "stop",
            "strand",
            "locus_tag",
            "gene",
            "product",
        ],
    )

    print(f"Wrote representative outputs to {args.out_dir.resolve()}")


if __name__ == "__main__":
    main()
