#!/usr/bin/env python3
"""Score all Yu PRSs directly from a Windows-mounted NAS without genotype copies."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_batched_prs_score as batch  # noqa: E402


def now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def run_command(command: list[str], label: str) -> None:
    print(f"{now()} | START {label}", flush=True)
    print("COMMAND | " + subprocess.list2cmdline(command), flush=True)
    completed = subprocess.run(command, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"{label} failed with exit code {completed.returncode}")
    print(f"{now()} | DONE {label}", flush=True)


def count_score_variants(score_file: Path, thresholds: list[tuple[str, float]]) -> list[int]:
    counts = [0] * len(thresholds)
    with score_file.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            for index, (column, _) in enumerate(thresholds):
                if float(row[column]) != 0:
                    counts[index] += 1
    return counts


def clump_ids(path: Path) -> set[str]:
    if not path.is_file() or path.stat().st_size == 0:
        return set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "ID" not in reader.fieldnames:
            raise ValueError(f"Clump output lacks ID column: {path}")
        return {row["ID"] for row in reader if row.get("ID")}


def write_wide_score(
    gwas: Path,
    ids: set[str],
    thresholds: list[tuple[str, float]],
    output: Path,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    seen: set[str] = set()
    with gzip.open(gwas, "rt", encoding="utf-8", newline="") as source, temporary.open(
        "w", encoding="utf-8", newline=""
    ) as target:
        reader = csv.DictReader(source, delimiter="\t")
        required = {"ID", "A1", "P", "BETA"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Normalized GWAS header is invalid: {gwas}")
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        writer.writerow(["ID", "A1", *[column for column, _ in thresholds]])
        for row in reader:
            variant_id = row["ID"]
            if variant_id not in ids:
                continue
            if variant_id in seen:
                raise ValueError(f"Duplicated clumped GWAS ID {variant_id}: {gwas}")
            seen.add(variant_id)
            p_value = float(row["P"])
            beta = float(row["BETA"])
            writer.writerow(
                [variant_id, row["A1"], *[beta if p_value <= threshold else 0 for _, threshold in thresholds]]
            )
    os.replace(temporary, output)


class DirectNasScorer:
    def __init__(self, args: argparse.Namespace):
        self.args = args
        self.project = Path(args.project_dir)
        self.analysis = Path(args.analysis_dir)
        self.prs = self.analysis / "13_prs"
        self.gwas = self.prs / "gwas_normalized"
        self.scores = self.prs / "scores"
        self.stream = self.prs / "genotype_stream"
        self.keep = self.prs / "inputs" / "target.keep.tsv"
        self.source_file = self.project / "config" / "prs_gwas_sources.tsv"
        self.threshold_file = self.project / "config" / "prs_thresholds.tsv"
        self.helper = self.project / "tools" / "build_batched_prs_score.py"
        self.sources = read_tsv(self.source_file)
        self.thresholds = batch.read_thresholds(self.threshold_file)
        self.nas = Path(args.nas_root)
        self.plink2 = str(Path(args.plink2))
        self.job_threads = max(1, args.workers // args.score_jobs)
        self.job_memory = max(2048, args.memory_mb // args.score_jobs)

    def validate(self) -> None:
        if len(self.sources) != 13 or len(self.thresholds) != 5:
            raise ValueError(f"Expected 13 outcomes and 5 thresholds; got {len(self.sources)} and {len(self.thresholds)}")
        for path in (Path(self.plink2), self.keep, self.helper):
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(f"Missing required file: {path}")
        for chromosome in range(self.args.start_chr, self.args.end_chr + 1):
            prefix = self.genotype_prefix(chromosome)
            for extension in ("pgen", "pvar", "psam"):
                path = Path(f"{prefix}.{extension}")
                if not path.is_file() or path.stat().st_size == 0:
                    raise FileNotFoundError(f"Missing direct NAS genotype: {path}")
            for source in self.sources:
                gwas = self.gwas / source["outcome_id"] / f"chr{chromosome}.tsv.gz"
                if not gwas.is_file() or gwas.stat().st_size == 0:
                    raise FileNotFoundError(f"Missing normalized GWAS: {gwas}")

    def genotype_prefix(self, chromosome: int) -> Path:
        return self.nas / f"chr{chromosome}" / "pgen" / f"chr{chromosome}_imp"

    def ensure_unique_ids(self, chromosome: int, genotype: Path, batch_dir: Path) -> Path:
        prefix = batch_dir / f"chr{chromosome}.unique"
        output = Path(f"{prefix}.snplist")
        if self.args.resume and output.is_file() and output.stat().st_size > 0:
            print(f"{now()} | REUSE unique IDs chr{chromosome}", flush=True)
            return output
        run_command(
            [
                self.plink2, "--pfile", str(genotype), "--rm-dup", "exclude-all",
                "--write-snplist", "--threads", str(self.job_threads),
                "--memory", str(self.job_memory), "--out", str(prefix),
            ],
            f"unique IDs chr{chromosome}",
        )
        if not output.is_file() or output.stat().st_size == 0:
            raise RuntimeError(f"Unique-ID list was not produced: {output}")
        return output

    def process_chromosome(self, chromosome: int) -> tuple[Path, Path]:
        genotype = self.genotype_prefix(chromosome)
        batch_dir = self.scores / "_batched" / f"chr{chromosome}"
        batch_dir.mkdir(parents=True, exist_ok=True)
        unique_ids = self.ensure_unique_ids(chromosome, genotype, batch_dir)
        status_rows: list[dict[str, object]] = []
        variant_rows: list[dict[str, object]] = []
        manifest_rows: list[dict[str, str]] = []

        for source in self.sources:
            outcome = source["outcome_id"]
            outdir = self.scores / outcome / f"chr{chromosome}"
            outdir.mkdir(parents=True, exist_ok=True)
            prefix = outdir / f"{outcome}_chr{chromosome}"
            sscore = Path(f"{prefix}.sscore")
            marker = Path(f"{prefix}.done")
            score_file = Path(f"{prefix}.wide_score.tsv")
            if self.args.resume and sscore.is_file() and marker.is_file():
                counts = count_score_variants(score_file, self.thresholds)
                for (column, threshold), count in zip(self.thresholds, counts):
                    variant_rows.append({
                        "outcome_id": outcome, "chr": chromosome, "threshold": threshold,
                        "score_column": column, "variant_count": count,
                    })
                status_rows.append({
                    "outcome_id": outcome, "chr": chromosome, "status": "PASS",
                    "sscore": str(sscore), "completed_at": now(),
                })
                print(f"{now()} | REUSE completed {outcome} chr{chromosome}", flush=True)
                continue

            gwas = self.gwas / outcome / f"chr{chromosome}.tsv.gz"
            clump_prefix = Path(f"{prefix}.clump")
            clumps = Path(f"{clump_prefix}.clumps")
            ids_file = Path(f"{prefix}.clumped.ids")
            if not (self.args.resume and score_file.is_file() and ids_file.is_file()):
                run_command(
                    [
                        self.plink2, "--pfile", str(genotype), "--extract", str(unique_ids),
                        "--clump", str(gwas), "--clump-id-field", "ID", "--clump-p-field", "P",
                        "--clump-a1-field", "A1", "--clump-force-a1",
                        "--clump-p1", str(self.args.max_p), "--clump-p2", str(self.args.max_p),
                        "--clump-r2", str(self.args.clump_r2), "--clump-kb", str(self.args.clump_kb),
                        "--threads", str(self.job_threads), "--memory", str(self.job_memory),
                        "--out", str(clump_prefix),
                    ],
                    f"clump {outcome} chr{chromosome}",
                )
                ids = clump_ids(clumps)
                ids_file.write_text("".join(f"{variant}\n" for variant in sorted(ids)), encoding="utf-8")
                write_wide_score(gwas, ids, self.thresholds, score_file)
            else:
                print(f"{now()} | REUSE clumping/weights {outcome} chr{chromosome}", flush=True)
            manifest_rows.append({
                "outcome_id": outcome, "score_file": str(score_file),
                "output_sscore": str(sscore), "marker": str(marker),
            })

        manifest = batch_dir / "batch_manifest.tsv"
        write_tsv(manifest, ["outcome_id", "score_file", "output_sscore", "marker"], manifest_rows)
        if manifest_rows:
            combined_score = batch_dir / f"chr{chromosome}.combined_score.tsv"
            batch_qc = batch_dir / "batch_qc.tsv"
            unresolved = batch_dir / "unresolved_alleles.tsv"
            batch.build(SimpleNamespace(
                thresholds=str(self.threshold_file), manifest=str(manifest),
                pvar=f"{genotype}.pvar", output=str(combined_score),
                qc=str(batch_qc), unresolved=str(unresolved),
            ))
            combined_prefix = batch_dir / f"chr{chromosome}.combined"
            combined_sscore = Path(f"{combined_prefix}.sscore")
            last_column = 2 + len(manifest_rows) * len(self.thresholds)
            run_command(
                [
                    self.plink2, "--pfile", str(genotype), "--keep", str(self.keep),
                    "--score", str(combined_score), "1", "2", "header-read", "ignore-dup-ids",
                    "cols=scoresums", "--score-col-nums", f"3-{last_column}",
                    "--threads", str(self.job_threads), "--memory", str(self.job_memory),
                    "--out", str(combined_prefix),
                ],
                f"batch score chr{chromosome}",
            )
            batch.split(SimpleNamespace(
                thresholds=str(self.threshold_file), manifest=str(manifest), input=str(combined_sscore)
            ))
            qc_rows = read_tsv(batch_qc)
            for row in qc_rows:
                variant_rows.append({
                    "outcome_id": row["outcome_id"], "chr": chromosome,
                    "threshold": row["threshold"], "score_column": row["score_column"],
                    "variant_count": row["output_nonzero"],
                })
            for row in manifest_rows:
                output = Path(row["output_sscore"])
                if not output.is_file() or output.stat().st_size == 0:
                    raise RuntimeError(f"Split score is missing: {output}")
                Path(row["marker"]).write_text(f"PASS\t{now()}\n", encoding="utf-8")
                status_rows.append({
                    "outcome_id": row["outcome_id"], "chr": chromosome, "status": "PASS",
                    "sscore": str(output), "completed_at": now(),
                })
            (batch_dir / "batch.done").write_text(
                f"PASS\t{now()}\toutcomes={len(manifest_rows)}\n", encoding="utf-8"
            )
            combined_sscore.unlink(missing_ok=True)

        status_path = batch_dir / "score_status.tsv"
        variants_path = batch_dir / "score_variant_counts.tsv"
        write_tsv(status_path, ["outcome_id", "chr", "status", "sscore", "completed_at"], status_rows)
        write_tsv(
            variants_path,
            ["outcome_id", "chr", "threshold", "score_column", "variant_count"],
            variant_rows,
        )
        print(f"{now()} | CHROMOSOME_PASS chr{chromosome} outcomes={len(status_rows)}", flush=True)
        return status_path, variants_path

    def refresh_global(self, results: list[tuple[Path, Path]]) -> None:
        status_path = self.prs / "score_status.tsv"
        variants_path = self.prs / "score_variant_counts.tsv"
        handled = set(range(self.args.start_chr, self.args.end_chr + 1))
        old_status = read_tsv(status_path) if status_path.is_file() else []
        old_variants = read_tsv(variants_path) if variants_path.is_file() else []
        status_rows = [row for row in old_status if int(row["chr"]) not in handled]
        variant_rows = [row for row in old_variants if int(row["chr"]) not in handled]
        for status_file, variant_file in results:
            status_rows.extend(read_tsv(status_file))
            variant_rows.extend(read_tsv(variant_file))
        order = {row["outcome_id"]: index for index, row in enumerate(self.sources)}
        status_rows.sort(key=lambda row: (int(row["chr"]), order[row["outcome_id"]]))
        variant_rows.sort(key=lambda row: (
            int(row["chr"]), order[row["outcome_id"]], float(row["threshold"])
        ))
        write_tsv(status_path, ["outcome_id", "chr", "status", "sscore", "completed_at"], status_rows)
        write_tsv(
            variants_path,
            ["outcome_id", "chr", "threshold", "score_column", "variant_count"],
            variant_rows,
        )

    def run(self) -> None:
        self.validate()
        print(
            f"{now()} | DIRECT_NAS_PASS root={self.nas} chromosomes={self.args.start_chr}-{self.args.end_chr} "
            f"jobs={self.args.score_jobs} threads/job={self.job_threads} memory/job={self.job_memory}",
            flush=True,
        )
        results: list[tuple[Path, Path]] = []
        with ThreadPoolExecutor(max_workers=self.args.score_jobs) as executor:
            futures = {
                executor.submit(self.process_chromosome, chromosome): chromosome
                for chromosome in range(self.args.start_chr, self.args.end_chr + 1)
            }
            for future in as_completed(futures):
                chromosome = futures[future]
                try:
                    results.append(future.result())
                except Exception as error:
                    raise RuntimeError(f"Chromosome {chromosome} failed: {error}") from error
        self.refresh_global(results)
        print(f"{now()} | PRS_SCORE_PASS chromosomes={len(results)}", flush=True)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument("--project-dir", required=True)
    result.add_argument("--analysis-dir", required=True)
    result.add_argument("--nas-root", required=True)
    result.add_argument("--plink2", required=True)
    result.add_argument("--workers", type=int, default=16)
    result.add_argument("--score-jobs", type=int, default=2)
    result.add_argument("--memory-mb", type=int, default=48000)
    result.add_argument("--start-chr", type=int, default=1)
    result.add_argument("--end-chr", type=int, default=22)
    result.add_argument("--max-p", type=float, default=0.0005)
    result.add_argument("--clump-r2", type=float, default=0.1)
    result.add_argument("--clump-kb", type=int, default=250)
    result.add_argument("--resume", action="store_true")
    return result


def main() -> int:
    args = parser().parse_args()
    if args.score_jobs < 1 or args.workers < 1 or args.memory_mb < 2048:
        print("ERROR: invalid resource parameters", file=sys.stderr)
        return 2
    try:
        DirectNasScorer(args).run()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
