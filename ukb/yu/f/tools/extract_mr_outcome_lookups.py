#!/usr/bin/env python3
"""Extract only local pQTL instrument variants from the 13 full outcome GWAS files."""

import argparse
import csv
import gzip
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--union-variants", required=True)
    p.add_argument("--gwas-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--workers", type=int, default=4)
    return p.parse_args()


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_targets(path):
    rsids, keys = set(), set()
    with open(path, newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            snp = row.get("snp", "").split(",", 1)[0].strip()
            chrom = row.get("chr", "").replace("chr", "").strip()
            pos = row.get("pos", "").split(".", 1)[0].strip()
            if snp and snp != ".":
                rsids.add(snp)
            if chrom and pos:
                keys.add(f"{chrom}:{pos}")
    return rsids, keys


def endpoint_from_name(path):
    return path.name.split(".", 1)[0]


def first_value(row, names, default=""):
    for name in names:
        value = row.get(name)
        if value not in (None, "", "NA", "nan"):
            return value
    return default


def normalize_row(row, outcome_id):
    if "#chrom" in row:  # FinnGen
        chrom = row["#chrom"].replace("chr", "")
        pos = row["pos"]
        snp = row.get("rsids", "").split(",", 1)[0]
        return {
            "outcome_id": outcome_id, "snp": snp, "chr": chrom, "pos": pos,
            "effect_allele_outcome": row["alt"].upper(),
            "other_allele_outcome": row["ref"].upper(),
            "eaf_outcome": row.get("af_alt", ""), "beta_outcome": row["beta"],
            "se_outcome": row["sebeta"], "p_outcome": row["pval"],
            "source_format": "FinnGen"
        }
    chrom = first_value(row, ["hm_chrom", "chromosome"]).replace("chr", "")
    pos = first_value(row, ["hm_pos", "base_pair_location"])
    return {
        "outcome_id": outcome_id,
        "snp": first_value(row, ["hm_rsid", "variant_id"]), "chr": chrom, "pos": pos,
        "effect_allele_outcome": first_value(row, ["hm_effect_allele", "effect_allele"]).upper(),
        "other_allele_outcome": first_value(row, ["hm_other_allele", "other_allele"]).upper(),
        "eaf_outcome": first_value(row, ["hm_effect_allele_frequency", "effect_allele_frequency"]),
        "beta_outcome": first_value(row, ["hm_beta", "beta"]),
        "se_outcome": first_value(row, ["standard_error"]),
        "p_outcome": first_value(row, ["p_value"]), "source_format": "GWAS_Catalog"
    }


def extract_one(path, rsids, keys, out_dir):
    endpoint = endpoint_from_name(path)
    out_path = out_dir / f"{endpoint}.csv.gz"
    fields = [
        "outcome_id", "snp", "chr", "pos", "effect_allele_outcome",
        "other_allele_outcome", "eaf_outcome", "beta_outcome", "se_outcome",
        "p_outcome", "source_format"
    ]
    kept = 0
    with gzip.open(path, "rt", newline="", encoding="utf-8-sig") as src, \
            gzip.open(out_path, "wt", newline="", encoding="utf-8") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        writer = csv.DictWriter(dst, fieldnames=fields)
        writer.writeheader()
        for raw in reader:
            row = normalize_row(raw, endpoint)
            snp = row["snp"].split(",", 1)[0]
            key = f"{row['chr']}:{str(row['pos']).split('.', 1)[0]}"
            if snp in rsids or key in keys:
                row["snp"] = snp
                writer.writerow(row)
                kept += 1
    return {
        "outcome_id": endpoint, "source_file": str(path), "source_sha256": sha256(path),
        "output_file": str(out_path), "rows": kept
    }


def main():
    args = parse_args()
    union = Path(args.union_variants)
    gwas_dir = Path(args.gwas_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rsids, keys = load_targets(union)
    files = sorted(gwas_dir.glob("*.gz"))
    if len(files) != 13:
        raise SystemExit(f"Expected 13 full outcome GWAS files, found {len(files)} in {gwas_dir}")
    results = []
    with ThreadPoolExecutor(max_workers=max(1, min(args.workers, 6))) as pool:
        futures = {pool.submit(extract_one, f, rsids, keys, out_dir): f for f in files}
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(f"{result['outcome_id']}: {result['rows']} rows", flush=True)
    results.sort(key=lambda x: x["outcome_id"])
    manifest = {
        "status": "PASS", "union_variants": str(union), "union_sha256": sha256(union),
        "target_rsids": len(rsids), "target_chr_pos": len(keys), "outcomes": results
    }
    with open(out_dir / "outcome_lookup_manifest.json", "w", encoding="utf-8") as f:
        json.dump(manifest, f, ensure_ascii=False, indent=2)
    if any(x["rows"] == 0 for x in results):
        raise SystemExit("One or more outcome lookups have zero matched variants; inspect build/ID harmonisation.")


if __name__ == "__main__":
    main()
