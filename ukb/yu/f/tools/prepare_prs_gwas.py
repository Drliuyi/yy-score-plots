#!/usr/bin/env python3
"""Download, normalize and split the 13 article-cited PRS GWAS files."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import re
import shutil
import urllib.request
from datetime import datetime
from pathlib import Path


def now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def valid_allele(value: str) -> bool:
    return bool(re.fullmatch(r"[ACGT]+", value.upper()))


def download(url: str, target: Path, resume: bool) -> None:
    partial = target.with_suffix(target.suffix + ".part")
    offset = partial.stat().st_size if resume and partial.exists() else 0
    request = urllib.request.Request(url, headers={"Range": f"bytes={offset}-"} if offset else {})
    print(f"{now()} | DOWNLOAD {url}", flush=True)
    with urllib.request.urlopen(request) as response, partial.open("ab" if offset and response.status == 206 else "wb") as output:
        shutil.copyfileobj(response, output, length=1024 * 1024)
    os.replace(partial, target)


def normalize(raw: Path, output: Path, source_kind: str, max_p: float) -> int:
    best: dict[str, tuple[int, int, str, str, float, float]] = {}
    with gzip.open(raw, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                if source_kind == "finngen_r9":
                    chromosome = int(row["#chrom"]); position = int(row["pos"])
                    ref, alt = row["ref"].upper(), row["alt"].upper()
                    variant = re.split(r"[,;]", row["rsids"])[0]
                    p_value, beta = float(row["pval"]), float(row["beta"])
                    a1, a2 = alt, ref
                elif source_kind == "gwas_catalog":
                    chromosome = int(row["hm_chrom"]); position = int(row["hm_pos"])
                    variant = row["hm_rsid"]
                    a1, a2 = row["hm_effect_allele"].upper(), row["hm_other_allele"].upper()
                    p_value, beta = float(row["p_value"]), float(row["hm_beta"])
                else:
                    raise ValueError(f"Unsupported GWAS source: {source_kind}")
            except (KeyError, TypeError, ValueError):
                continue
            if not (1 <= chromosome <= 22 and re.fullmatch(r"rs\d+", variant or "") and 0 < p_value <= max_p and valid_allele(a1) and valid_allele(a2)):
                continue
            candidate = (chromosome, position, a1, a2, beta, p_value)
            if variant not in best or p_value < best[variant][-1]:
                best[variant] = candidate
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with gzip.open(temporary, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["ID", "CHR", "POS", "A1", "A2", "BETA", "P"])
        for variant, values in sorted(best.items(), key=lambda item: (item[1][0], item[0], item[1][-1])):
            writer.writerow([variant, *values])
    os.replace(temporary, output)
    return len(best)


def split_chromosomes(normalized: Path, output_dir: Path, resume: bool) -> None:
    handles = {}
    writers = {}
    try:
        for chromosome in range(1, 23):
            target = output_dir / f"chr{chromosome}.tsv.gz"
            if resume and target.is_file() and target.stat().st_size:
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            temporary = target.with_suffix(target.suffix + ".tmp")
            handles[chromosome] = gzip.open(temporary, "wt", encoding="utf-8", newline="")
            writers[chromosome] = csv.writer(handles[chromosome], delimiter="\t", lineterminator="\n")
            writers[chromosome].writerow(["ID", "CHR", "POS", "A1", "A2", "BETA", "P"])
        with gzip.open(normalized, "rt", encoding="utf-8", newline="") as source:
            for row in csv.DictReader(source, delimiter="\t"):
                chromosome = int(row["CHR"])
                if chromosome in writers:
                    writers[chromosome].writerow([row[key] for key in ("ID", "CHR", "POS", "A1", "A2", "BETA", "P")])
    finally:
        for handle in handles.values():
            handle.close()
    for chromosome in handles:
        target = output_dir / f"chr{chromosome}.tsv.gz"
        os.replace(target.with_suffix(target.suffix + ".tmp"), target)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", required=True)
    parser.add_argument("--analysis-dir", required=True)
    parser.add_argument("--max-p", type=float, default=0.0005)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    project = Path(args.project_dir)
    prs = Path(args.analysis_dir) / "13_prs"
    raw_dir, normalized_dir = prs / "gwas_raw", prs / "gwas_normalized"
    raw_dir.mkdir(parents=True, exist_ok=True); normalized_dir.mkdir(parents=True, exist_ok=True)
    with (project / "config/prs_gwas_sources.tsv").open(encoding="utf-8-sig", newline="") as handle:
        sources = list(csv.DictReader(handle, delimiter="\t"))
    statuses = []
    for source in sources:
        outcome, source_id = source["outcome_id"], source["source_id"]
        raw = raw_dir / f"{outcome}.{source_id}.gz"
        normalized = normalized_dir / f"{outcome}.maxp_{args.max_p}.tsv.gz"
        print(f"{now()} | {outcome}", flush=True)
        if not (args.resume and normalized.is_file() and normalized.stat().st_size):
            if not raw.is_file() or not raw.stat().st_size:
                download(source["url"], raw, args.resume)
            retained = normalize(raw, normalized, source["source_kind"], args.max_p)
        else:
            with gzip.open(normalized, "rt", encoding="utf-8") as handle:
                retained = sum(1 for _ in handle) - 1
        if retained < 1:
            raise RuntimeError(f"No variants retained for {outcome}")
        split_chromosomes(normalized, normalized_dir / outcome, args.resume)
        statuses.append({
            "outcome_id": outcome, "source_id": source_id, "status": "PASS",
            "raw_file": str(raw), "normalized_file": str(normalized),
            "retained_variants": retained, "sha256": sha256(normalized), "completed_at": now(),
        })
    status = prs / "gwas_prepare_status.tsv"
    with status.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(statuses[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(statuses)
    print(f"GWAS PREPARATION PASS | outcomes={len(statuses)}")


if __name__ == "__main__":
    main()
