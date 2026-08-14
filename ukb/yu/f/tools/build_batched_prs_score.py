#!/usr/bin/env python3
"""Build and split one-pass, multi-outcome PLINK2 score files."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import sys
from collections import defaultdict
from pathlib import Path


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def norm_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]", "", value).upper()


def read_thresholds(path: Path) -> list[tuple[str, float]]:
    with open_text(path) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 5:
        raise ValueError(f"Expected five PRS thresholds, found {len(rows)}")
    return [(row["score_column"], float(row["threshold"])) for row in rows]


def read_manifest(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"outcome_id", "score_file", "output_sscore"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Batch manifest must contain {sorted(required)}")
    outcomes = [row["outcome_id"] for row in rows]
    if len(outcomes) != len(set(outcomes)):
        raise ValueError("Batch manifest has duplicated outcome_id values")
    return rows


def read_score_files(
    manifest: list[dict[str, str]], thresholds: list[tuple[str, float]]
):
    records: dict[str, dict[str, tuple[str, list[float]]]] = defaultdict(dict)
    source_counts: dict[tuple[str, str], int] = defaultdict(int)
    score_columns = [item[0] for item in thresholds]

    for row in manifest:
        outcome = row["outcome_id"]
        path = Path(row["score_file"])
        if not path.is_file():
            raise FileNotFoundError(f"Missing outcome score file: {path}")
        with open_text(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"ID", "A1", *score_columns}
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise ValueError(f"Score file has invalid header: {path}")
            for source in reader:
                variant_id = source["ID"].strip()
                if not variant_id or variant_id == ".":
                    continue
                if outcome in records[variant_id]:
                    raise ValueError(f"Duplicated {outcome}/{variant_id} in {path}")
                weights = [float(source[column]) for column in score_columns]
                records[variant_id][outcome] = (source["A1"].upper(), weights)
                for (column, _), weight in zip(thresholds, weights):
                    if weight != 0:
                        source_counts[(outcome, column)] += 1
    return records, source_counts


def read_needed_pvar_alleles(path: Path, needed: set[str]):
    alleles: dict[str, tuple[str, str]] = {}
    duplicate_ids: set[str] = set()
    header = None
    indexes = None

    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if header is None:
                header = [field.lstrip("#") for field in fields]
                required = {"ID", "REF", "ALT"}
                if not required.issubset(header):
                    raise ValueError(f"PVAR header is missing {sorted(required)}: {path}")
                indexes = {name: header.index(name) for name in required}
                continue
            assert indexes is not None
            variant_id = fields[indexes["ID"]]
            if variant_id not in needed:
                continue
            if variant_id in alleles:
                duplicate_ids.add(variant_id)
                continue
            alleles[variant_id] = (
                fields[indexes["REF"]].upper(), fields[indexes["ALT"]].upper()
            )

    for variant_id in duplicate_ids:
        alleles.pop(variant_id, None)
    return alleles, duplicate_ids


def build(args: argparse.Namespace) -> None:
    thresholds = read_thresholds(Path(args.thresholds))
    manifest = read_manifest(Path(args.manifest))
    records, source_counts = read_score_files(manifest, thresholds)
    alleles, duplicate_ids = read_needed_pvar_alleles(Path(args.pvar), set(records))

    outcomes = [row["outcome_id"] for row in manifest]
    columns = [(outcome, score_column) for outcome in outcomes for score_column, _ in thresholds]
    output_counts: dict[tuple[str, str], int] = defaultdict(int)
    orientation_counts: dict[tuple[str, str], int] = defaultdict(int)
    unresolved: list[tuple[str, str, str, str, str]] = []
    rows_out: list[tuple[str, str, list[float]]] = []

    for variant_id in sorted(records):
        target = alleles.get(variant_id)
        if target is None:
            reason = "duplicate_target_id" if variant_id in duplicate_ids else "missing_target_id"
            for outcome, (a1, _) in records[variant_id].items():
                unresolved.append((outcome, variant_id, a1, "", reason))
            continue
        ref, alt = target
        if "," in alt:
            for outcome, (a1, _) in records[variant_id].items():
                unresolved.append((outcome, variant_id, a1, f"{ref}/{alt}", "multiallelic"))
            continue

        values = {(outcome, column): 0.0 for outcome, column in columns}
        matched_any = False
        for outcome, (a1, weights) in records[variant_id].items():
            if a1 == alt:
                sign = 1.0
                orientation_counts[(outcome, "same_as_alt")] += 1
            elif a1 == ref:
                sign = -1.0
                orientation_counts[(outcome, "flipped_from_ref")] += 1
            else:
                unresolved.append((outcome, variant_id, a1, f"{ref}/{alt}", "allele_mismatch"))
                continue
            matched_any = True
            for (column, _), weight in zip(thresholds, weights):
                oriented = sign * weight
                values[(outcome, column)] = oriented
                if oriented != 0:
                    output_counts[(outcome, column)] += 1
        if matched_any:
            rows_out.append((variant_id, alt, [values[column] for column in columns]))

    qc_path = Path(args.qc)
    qc_path.parent.mkdir(parents=True, exist_ok=True)
    with qc_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["outcome_id", "threshold", "score_column", "source_nonzero", "output_nonzero",
             "same_as_alt", "flipped_from_ref", "unresolved"]
        )
        unresolved_counts = defaultdict(int)
        for outcome, *_ in unresolved:
            unresolved_counts[outcome] += 1
        for outcome in outcomes:
            for column, threshold in thresholds:
                writer.writerow(
                    [outcome, threshold, column, source_counts[(outcome, column)],
                     output_counts[(outcome, column)], orientation_counts[(outcome, "same_as_alt")],
                     orientation_counts[(outcome, "flipped_from_ref")], unresolved_counts[outcome]]
                )

    unresolved_path = Path(args.unresolved)
    unresolved_path.parent.mkdir(parents=True, exist_ok=True)
    with unresolved_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["outcome_id", "ID", "source_A1", "target_REF_ALT", "reason"])
        writer.writerows(unresolved)

    if unresolved:
        raise RuntimeError(
            f"Refusing batched scoring: {len(unresolved)} score rows could not be aligned; "
            f"see {unresolved_path}"
        )
    if not rows_out:
        raise RuntimeError("No aligned variants remain for batched scoring")

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = output_path.with_suffix(output_path.suffix + ".tmp")
    with temp_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["ID", "A1", *[f"{outcome}__{column}" for outcome, column in columns]])
        for variant_id, alt, values in rows_out:
            writer.writerow([variant_id, alt, *[format(value, ".12g") for value in values]])
    os.replace(temp_path, output_path)
    print(
        f"BATCH_BUILD_PASS outcomes={len(outcomes)} columns={len(columns)} "
        f"variants={len(rows_out)} output={output_path}"
    )


def resolve_column(fieldnames: list[str], requested: str) -> str:
    target = norm_name(requested)
    hits = [field for field in fieldnames if norm_name(field) == target]
    if len(hits) != 1:
        raise ValueError(f"Cannot resolve exactly one column {requested}; hits={hits}")
    return hits[0]


def split(args: argparse.Namespace) -> None:
    thresholds = read_thresholds(Path(args.thresholds))
    manifest = read_manifest(Path(args.manifest))
    input_path = Path(args.input)
    outputs = []

    with open_text(input_path) as source:
        reader = csv.DictReader(source, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Empty PLINK score output: {input_path}")
        iid = resolve_column(reader.fieldnames, "IID")
        fid_hits = [field for field in reader.fieldnames if norm_name(field) == "FID"]
        if len(fid_hits) > 1:
            raise ValueError(f"Cannot resolve at most one FID column; hits={fid_hits}")
        fid = fid_hits[0] if fid_hits else None
        for row in manifest:
            outcome = row["outcome_id"]
            output_path = Path(row["output_sscore"])
            output_path.parent.mkdir(parents=True, exist_ok=True)
            temp_path = output_path.with_suffix(output_path.suffix + ".tmp")
            handle = temp_path.open("w", encoding="utf-8", newline="")
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            requested = [f"{outcome}__{column}_SUM" for column, _ in thresholds]
            observed = [resolve_column(reader.fieldnames, name) for name in requested]
            writer.writerow(["#FID", "IID", *[f"{column}_SUM" for column, _ in thresholds]])
            outputs.append((output_path, temp_path, handle, writer, observed))

        row_count = 0
        try:
            for source_row in reader:
                row_count += 1
                for _, _, _, writer, observed in outputs:
                    writer.writerow([
                        source_row[fid] if fid else source_row[iid], source_row[iid],
                        *[source_row[name] for name in observed]
                    ])
        finally:
            for _, _, handle, _, _ in outputs:
                handle.close()

    if row_count == 0:
        for _, temp_path, _, _, _ in outputs:
            temp_path.unlink(missing_ok=True)
        raise RuntimeError(f"PLINK score output has no participant rows: {input_path}")
    for output_path, temp_path, _, _, _ in outputs:
        os.replace(temp_path, output_path)
    print(f"BATCH_SPLIT_PASS outcomes={len(outputs)} participants={row_count}")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    subparsers = result.add_subparsers(dest="command", required=True)

    build_parser = subparsers.add_parser("build")
    build_parser.add_argument("--manifest", required=True)
    build_parser.add_argument("--thresholds", required=True)
    build_parser.add_argument("--pvar", required=True)
    build_parser.add_argument("--output", required=True)
    build_parser.add_argument("--qc", required=True)
    build_parser.add_argument("--unresolved", required=True)
    build_parser.set_defaults(func=build)

    split_parser = subparsers.add_parser("split")
    split_parser.add_argument("--manifest", required=True)
    split_parser.add_argument("--thresholds", required=True)
    split_parser.add_argument("--input", required=True)
    split_parser.set_defaults(func=split)
    return result


def main() -> int:
    args = parser().parse_args()
    try:
        args.func(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
