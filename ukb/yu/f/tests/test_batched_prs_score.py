#!/usr/bin/env python3
import csv
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]
HELPER = PROJECT_DIR / "tools" / "build_batched_prs_score.py"


def write_tsv(path: Path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)


class BatchedPrsScoreTest(unittest.TestCase):
    def test_build_orients_alleles_and_split_restores_contract(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            thresholds = root / "thresholds.tsv"
            write_tsv(
                thresholds,
                [
                    ["threshold", "label", "score_column"],
                    ["0.000005", "5e-06", "PT_5e_06"],
                    ["0.00001", "1e-05", "PT_1e_05"],
                    ["0.00005", "5e-05", "PT_5e_05"],
                    ["0.0001", "0.0001", "PT_1e_04"],
                    ["0.0005", "0.0005", "PT_5e_04"],
                ],
            )
            score_header = [
                "ID", "A1", "PT_5e_06", "PT_1e_05", "PT_5e_05", "PT_1e_04", "PT_5e_04"
            ]
            score_a = root / "a.tsv"
            score_b = root / "b.tsv"
            write_tsv(score_a, [score_header, ["rs1", "A", 1, 2, 3, 4, 5], ["rs2", "C", 0, 1, 2, 3, 4]])
            write_tsv(score_b, [score_header, ["rs1", "G", 5, 4, 3, 2, 1], ["rs3", "T", 1, 1, 1, 1, 1]])

            output_a = root / "a.sscore"
            output_b = root / "b.sscore"
            manifest = root / "manifest.tsv"
            write_tsv(
                manifest,
                [
                    ["outcome_id", "score_file", "output_sscore", "marker"],
                    ["outcome_a", score_a, output_a, root / "a.done"],
                    ["outcome_b", score_b, output_b, root / "b.done"],
                ],
            )
            pvar = root / "chr1.pvar"
            write_tsv(
                pvar,
                [
                    ["#CHROM", "POS", "ID", "REF", "ALT"],
                    [1, 100, "rs1", "G", "A"],
                    [1, 200, "rs2", "C", "T"],
                    [1, 300, "rs3", "C", "T"],
                ],
            )
            combined = root / "combined.tsv"
            qc = root / "qc.tsv"
            unresolved = root / "unresolved.tsv"
            subprocess.run(
                [
                    sys.executable, str(HELPER), "build", "--manifest", str(manifest),
                    "--thresholds", str(thresholds), "--pvar", str(pvar),
                    "--output", str(combined), "--qc", str(qc), "--unresolved", str(unresolved),
                ],
                check=True,
            )
            with combined.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            by_id = {row["ID"]: row for row in rows}
            self.assertEqual(float(by_id["rs1"]["outcome_a__PT_5e_06"]), 1.0)
            self.assertEqual(float(by_id["rs1"]["outcome_b__PT_5e_06"]), -5.0)
            self.assertEqual(float(by_id["rs2"]["outcome_a__PT_1e_05"]), -1.0)
            self.assertEqual(float(by_id["rs3"]["outcome_b__PT_5e_04"]), 1.0)

            combined_sscore = root / "combined.sscore"
            score_names = [field for field in by_id["rs1"] if field not in {"ID", "A1"}]
            plink2 = shutil.which("plink2")
            if plink2:
                vcf = root / "tiny.vcf"
                samples = [str(1000 + index) for index in range(1, 51)]
                genotypes = [["0/0", "0/1", "1/1"][index % 3] for index in range(50)]
                vcf.write_text("\n".join([
                    "##fileformat=VCFv4.2",
                    "##contig=<ID=1,length=1000000>",
                    "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples]),
                    "\t".join(["1", "100", "rs1", "G", "A", ".", "PASS", ".", "GT", *genotypes]),
                    "\t".join(["1", "200", "rs2", "C", "T", ".", "PASS", ".", "GT", *genotypes[1:] + genotypes[:1]]),
                    "\t".join(["1", "300", "rs3", "C", "T", ".", "PASS", ".", "GT", *genotypes[2:] + genotypes[:2]]),
                ]) + "\n", encoding="utf-8")
                geno = root / "tiny"
                subprocess.run([plink2, "--vcf", str(vcf), "--make-pgen", "--out", str(geno)], check=True)
                scored = root / "combined"
                subprocess.run(
                    [
                        plink2, "--pfile", str(geno), "--score", str(combined), "1", "2",
                        "header-read", "ignore-dup-ids", "cols=scoresums",
                        "--score-col-nums", f"3-{2 + len(score_names)}", "--out", str(scored),
                    ],
                    check=True,
                )
            else:
                write_tsv(
                    combined_sscore,
                    [
                        ["#FID", "IID", *[f"{field}_SUM" for field in score_names]],
                        ["1001", "1001", *range(1, len(score_names) + 1)],
                        ["1002", "1002", *range(101, 101 + len(score_names))],
                    ],
                )
            subprocess.run(
                [
                    sys.executable, str(HELPER), "split", "--manifest", str(manifest),
                    "--thresholds", str(thresholds), "--input", str(combined_sscore),
                ],
                check=True,
            )
            for output in (output_a, output_b):
                with output.open(encoding="utf-8") as handle:
                    split_rows = list(csv.reader(handle, delimiter="\t"))
                self.assertEqual(split_rows[0][0:2], ["#FID", "IID"])
                self.assertEqual(split_rows[0][2:], [
                    "PT_5e_06_SUM", "PT_1e_05_SUM", "PT_5e_05_SUM", "PT_1e_04_SUM", "PT_5e_04_SUM"
                ])
                self.assertEqual(len(split_rows), 51 if plink2 else 3)


if __name__ == "__main__":
    unittest.main()
