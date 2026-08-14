#!/usr/bin/env python3
"""Windows-only end-to-end smoke test for the direct-NAS PRS scorer."""

import csv
import gzip
import subprocess
import sys
import tempfile
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
SCORER = PROJECT / "tools" / "score_prs_directnas_windows.py"
SOURCES = PROJECT / "config" / "prs_gwas_sources.tsv"


def main() -> int:
    if sys.platform != "win32":
        print("SKIP: Windows-only integration test")
        return 0
    plink2 = PROJECT / "tools" / "bin" / "plink2.exe"
    if not plink2.is_file():
        raise FileNotFoundError(plink2)
    with SOURCES.open(encoding="utf-8-sig", newline="") as handle:
        outcomes = [row["outcome_id"] for row in csv.DictReader(handle, delimiter="\t")]
    if len(outcomes) != 13:
        raise AssertionError(f"Expected 13 outcomes, found {len(outcomes)}")

    with tempfile.TemporaryDirectory(dir=PROJECT / "tests") as temporary:
        root = Path(temporary)
        analysis = root / "analysis"
        nas = root / "nas"
        genotype_dir = nas / "chr1" / "pgen"
        genotype_dir.mkdir(parents=True)
        samples = [str(100000 + index) for index in range(60)]
        vcf = root / "chr1.vcf"
        lines = [
            "##fileformat=VCFv4.2",
            "##contig=<ID=1,length=1000000>",
            "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples]),
        ]
        genotypes = ["0/0", "0/1", "1/1"]
        for variant in range(1, 101):
            calls = [genotypes[(sample + variant * 7) % 3] for sample in range(len(samples))]
            lines.append("\t".join([
                "1", str(variant * 1000), f"rs{variant}", "G", "A", ".", "PASS", ".", "GT", *calls
            ]))
        vcf.write_text("\n".join(lines) + "\n", encoding="utf-8")
        genotype_prefix = genotype_dir / "chr1_imp"
        subprocess.run([
            str(plink2), "--vcf", str(vcf), "--double-id", "--make-pgen", "--out", str(genotype_prefix)
        ], check=True)

        inputs = analysis / "13_prs" / "inputs"
        inputs.mkdir(parents=True)
        with (inputs / "target.keep.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["#FID", "IID"])
            writer.writerows([[sample, sample] for sample in samples])
        for outcome_index, outcome in enumerate(outcomes):
            output_dir = analysis / "13_prs" / "gwas_normalized" / outcome
            output_dir.mkdir(parents=True)
            with gzip.open(output_dir / "chr1.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["ID", "A1", "P", "BETA"])
                for variant in range(1, 101):
                    writer.writerow([
                        f"rs{variant}", "A", variant * 4e-6,
                        (1 if (variant + outcome_index) % 2 else -1) * (0.01 + variant / 10000),
                    ])

        subprocess.run([
            sys.executable, str(SCORER), "--project-dir", str(PROJECT),
            "--analysis-dir", str(analysis), "--nas-root", str(nas),
            "--plink2", str(plink2), "--workers", "4", "--score-jobs", "1",
            "--memory-mb", "4096", "--start-chr", "1", "--end-chr", "1", "--resume",
        ], check=True)
        score_files = list((analysis / "13_prs" / "scores").glob("*/chr1/*.sscore"))
        if len(score_files) != 13:
            raise AssertionError(f"Expected 13 score files, found {len(score_files)}")
        with (analysis / "13_prs" / "score_status.tsv").open(encoding="utf-8") as handle:
            status_rows = list(csv.DictReader(handle, delimiter="\t"))
        with (analysis / "13_prs" / "score_variant_counts.tsv").open(encoding="utf-8") as handle:
            variant_rows = list(csv.DictReader(handle, delimiter="\t"))
        if len(status_rows) != 13 or len(variant_rows) != 65:
            raise AssertionError(f"Invalid output matrix: status={len(status_rows)} variants={len(variant_rows)}")
        print("DIRECTNAS_WINDOWS_E2E_PASS outcomes=13 thresholds=5 chromosomes=1")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
