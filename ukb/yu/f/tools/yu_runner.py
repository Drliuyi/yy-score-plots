#!/usr/bin/env python3
"""Single command runner for the frozen Yu/Chen proteomic workflow."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[2]
F_DIR = PROJECT / "f"
FULL_R = F_DIR / "entry" / "99_run_yu_full_reproduction.R"
PRS_R = F_DIR / "entry" / "99_run_yu_prs.R"
MODEL_PY = F_DIR / "python" / "04_full_reproduction.py"
MODEL_CONFIG = F_DIR / "config" / "full_reproduction_defaults.json"
GWAS_PREP = F_DIR / "tools" / "prepare_prs_gwas.py"
PRS_SCORER = F_DIR / "tools" / "score_prs_directnas_windows.py"
EXPECTED_XLSX = "e21e5699f2b9fb7f86d26fed90bbebe03688d46f4a035a1f6494cefc14b895d7"
EXPECTED_PDF = "b78b07054109605b15b99a0501b7527448b073e9ce8330addf7d637aa655f29a"

STEPS = {
    1: ("Sources and preflight", "light"),
    2: ("Incident cohort", "light"),
    3: ("Full-panel Cox associations", "heavy"),
    4: ("Protein selection and prediction", "heavy"),
    5: ("CMR associations", "heavy"),
    6: ("Mendelian randomization", "medium"),
    7: ("Mediation candidate analysis", "medium"),
    8: ("CMAverse mediation", "very heavy"),
    9: ("PRS-protein reconstruction", "very heavy"),
    10: ("Enrichment, TF and PPI", "medium/network"),
    11: ("Final Figures 1-6 and report", "light"),
}


def stamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def fail(message: str) -> None:
    raise RuntimeError(message)


def is_windows_path(value: str) -> bool:
    return bool(re.match(r"^[A-Za-z]:[\\/]", value or ""))


def local_path(value: str) -> str:
    value = str(value or "").strip().replace("\\", "/")
    if os.name == "nt" or not is_windows_path(value):
        return value.rstrip("/") if value not in ("/", "") else value
    if shutil.which("wslpath"):
        result = subprocess.run(["wslpath", "-u", value], text=True, capture_output=True)
        if result.returncode == 0:
            return result.stdout.strip()
        drive, rest = value[0].lower(), value[2:].lstrip("/")
        return f"/mnt/{drive}/{rest}"
    if shutil.which("cygpath"):
        return subprocess.check_output(["cygpath", "-u", value], text=True).strip()
    return value


def join_root(root: str, *parts: str) -> str:
    prefix = root.rstrip("/\\")
    if not prefix and root.startswith("/"):
        prefix = "/"
    suffix = "/".join(part.strip("/\\") for part in parts)
    return f"{prefix}/{suffix}" if prefix != "/" else f"/{suffix}"


def runtime_path(value: str | Path, executable: str) -> str:
    value = str(value)
    if not executable.lower().endswith(".exe") or os.name == "nt":
        return value
    if is_windows_path(value):
        return value.replace("/", "\\")
    converter = "wslpath" if shutil.which("wslpath") else "cygpath"
    if converter:
        return subprocess.check_output([converter, "-w", value], text=True).strip()
    return value


def existing_file(value: str) -> bool:
    return bool(value) and Path(local_path(value)).is_file()


def existing_dir(value: str) -> bool:
    return bool(value) and Path(local_path(value)).is_dir()


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def valid_model_runtime(output: str) -> bool:
    fields = output.strip().split()
    return len(fields) == 2 and fields[0].split(".")[:2] == ["3", "9"] and fields[1] == "3.3.2"


def profile_path(explicit: str) -> Path:
    if explicit:
        return Path(local_path(explicit)).expanduser()
    root = Path(os.environ.get("XDG_CONFIG_HOME", Path.home() / ".config"))
    return root / "yu-protein-analysis" / "paths.json"


def load_profile(path: Path) -> dict[str, str]:
    if not path.is_file():
        return {}
    try:
        content = json.loads(path.read_text(encoding="utf-8-sig"))
        return {key: local_path(str(value)) for key, value in content.items() if isinstance(value, str)}
    except Exception as error:
        fail(f"Invalid path profile: {path}. Delete it or use --reset-paths. {error}")


def save_profile(path: Path, values: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {"schema_version": 2, "updated_at": stamp(), **values}
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def first_existing(values: list[str], kind: str = "file") -> str:
    check = existing_file if kind == "file" else existing_dir
    for value in values:
        if value and check(value):
            return local_path(value)
    return ""


def win_user_roots() -> list[str]:
    roots: list[str] = []
    if os.name == "nt" and os.environ.get("USERPROFILE"):
        roots.append(os.environ["USERPROFILE"])
    users = Path("/mnt/c/Users")
    if users.is_dir():
        roots.extend(str(path) for path in users.iterdir() if path.is_dir())
    return roots


def resolve_executable(explicit: str, commands: list[str], known: list[str], label: str) -> str:
    candidates = [explicit] if explicit else []
    candidates.extend(known)
    for command in commands:
        found = shutil.which(command)
        if found:
            candidates.append(found)
    for candidate in candidates:
        path = local_path(candidate)
        if path and Path(path).is_file():
            return path
    fail(f"{label} not found. Run './yu.sh install' or supply its path override.")


def resolve_rscript(explicit: str) -> str:
    known = [
        "/mnt/c/Program Files/R/R-4.5.1/bin/x64/Rscript.exe",
        "/mnt/c/Program Files/R/R-4.3.2/bin/x64/Rscript.exe",
        "C:/Program Files/R/R-4.5.1/bin/x64/Rscript.exe",
        "C:/Program Files/R/R-4.3.2/bin/x64/Rscript.exe",
    ]
    return resolve_executable(explicit or os.environ.get("YU_RSCRIPT", ""), ["Rscript"], known, "Rscript")


def resolve_model_python(explicit: str) -> str:
    known: list[str] = []
    for root in win_user_roots():
        known.extend([
            f"{root}/anaconda3/envs/yu_proteomic_repo_py39/python.exe",
            f"{root}/miniconda3/envs/yu_proteomic_repo_py39/python.exe",
        ])
    known.extend([
        "/mnt/c/ProgramData/anaconda3/envs/yu_proteomic_repo_py39/python.exe",
        "/mnt/c/ProgramData/miniconda3/envs/yu_proteomic_repo_py39/python.exe",
    ])
    return resolve_executable(explicit or os.environ.get("YU_PYTHON", ""), [], known, "Frozen Python 3.9")


def resolve_conda(explicit: str) -> str:
    known = [
        "/mnt/c/ProgramData/anaconda3/Scripts/conda.exe",
        "/mnt/c/ProgramData/miniconda3/Scripts/conda.exe",
        "C:/ProgramData/anaconda3/Scripts/conda.exe",
        "C:/ProgramData/miniconda3/Scripts/conda.exe",
    ]
    for root in win_user_roots():
        known.extend([f"{root}/anaconda3/Scripts/conda.exe", f"{root}/miniconda3/Scripts/conda.exe"])
    return resolve_executable(explicit or os.environ.get("YU_CONDA", ""), ["conda"], known, "Conda")


def tee_run(command: list[str], log: Path, env: dict[str, str] | None = None) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    print(f"{stamp()} | LOG {log}", flush=True)
    with log.open("w", encoding="utf-8") as handle:
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, encoding="utf-8", errors="replace", env=env,
        )
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="", flush=True)
            handle.write(line)
        code = process.wait()
    if code:
        fail(f"Command failed with exit code {code}. See {log}")


def read_first_column(path: Path, delimiter: str = ",") -> list[str]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = csv.DictReader(handle, delimiter=delimiter)
        key = rows.fieldnames[0] if rows.fieldnames else ""
        return [row[key].strip() for row in rows if row.get(key, "").strip()]


def parse_steps(value: str) -> list[int]:
    aliases = {
        "core": list(range(1, 5)), "downstream": list(range(5, 11)),
        "all": list(range(1, 12)), "figures": [11], "finalize": [10, 11],
    }
    value = value.strip().lower()
    if value in aliases:
        return aliases[value]
    result: set[int] = set()
    for token in re.split(r"[,; ]+", value):
        if not token:
            continue
        match = re.fullmatch(r"(\d+)-(\d+)", token)
        if match:
            start, end = map(int, match.groups())
            if start > end:
                fail(f"Invalid descending step range: {token}")
            result.update(range(start, end + 1))
        elif token.isdigit():
            result.add(int(token))
        else:
            fail(f"Unknown step: {token}")
    if not result or min(result) < 1 or max(result) > 11:
        fail("Steps must be between 1 and 11.")
    return sorted(result)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="./yu.sh",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Yu Protein Analysis

Commands:
  ./yu.sh setup                 save local paths; no analysis
  ./yu.sh doctor                check code, data and runtimes
  ./yu.sh 1-4 --resume          cohort, Cox and prediction
  ./yu.sh all --confirm-heavy   complete Steps 1-11
  ./yu.sh finalize --resume     systems analysis and final figures
  ./yu.sh status                read-only progress
  ./yu.sh package               source ZIP and SHA256

Steps:
  1 sources/QC     2 cohort       3 Cox          4 prediction
  5 CMR            6 MR           7 mediation    8 CMAverse
  9 PRS            10 systems     11 figures/report

Steps 8 and 9 require --confirm-heavy. Disease-specific or custom-panel runs
require a new project.
""",
    )
    parser.add_argument("command", nargs="?", default="help")
    parser.add_argument("--step", dest="step")
    parser.add_argument("--analysis-project", default="yu_proteomic_repo_v3")
    parser.add_argument("--disease", default="all")
    parser.add_argument("--protein-panel", choices=["local_reselected", "published_257", "custom"], default="local_reselected")
    parser.add_argument("--model-protein-file", default="")
    parser.add_argument("--model-proteins", default="")
    parser.add_argument("--figure4-extra-project", default="")
    parser.add_argument("--figure4-extra-outcome", default="")
    parser.add_argument("--figure4-extra-label", default="")
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--cox-jobs", type=int, default=4)
    parser.add_argument("--cmr-jobs", type=int, default=4)
    parser.add_argument("--model-jobs", type=int, default=3)
    parser.add_argument("--cmest-jobs", type=int, default=8)
    parser.add_argument("--association-jobs", type=int, default=4)
    parser.add_argument("--score-jobs", type=int, default=2)
    parser.add_argument("--bootstrap-n", type=int, default=1000)
    parser.add_argument("--cmest-pilot-boot", type=int, default=20)
    parser.add_argument("--memory-mb", type=int, default=48000)
    parser.add_argument("--dir0", default=os.environ.get("YU_DIR0", ""))
    parser.add_argument("--phe-dir", default=os.environ.get("YU_PHEDIR", ""))
    parser.add_argument("--analysis-root", default=os.environ.get("YU_OUTDIR", ""))
    parser.add_argument("--raw-protein-file", default="")
    parser.add_argument("--phenotype-rds", default="")
    parser.add_argument("--raw-phenotype-file", default="")
    parser.add_argument("--panel-mapping-file", default="")
    parser.add_argument("--cmr-feature-file", default="")
    parser.add_argument("--pqtl-root", default="")
    parser.add_argument("--mr-outcome-lookup-dir", default="")
    parser.add_argument("--supplement-workbook-file", default="")
    parser.add_argument("--supplement-methods-file", default="")
    parser.add_argument("--olink-processing-start-date-file", default="")
    parser.add_argument("--genotype-root", default="")
    parser.add_argument("--rscript", default="")
    parser.add_argument("--model-python", default="")
    parser.add_argument("--conda", default=os.environ.get("YU_CONDA", ""))
    parser.add_argument("--path-config", default=os.environ.get("YU_PATH_CONFIG", ""))
    parser.add_argument("--path-prompt", choices=["auto", "console", "off"], default="auto")
    parser.add_argument("--systems-score", type=int, default=700)
    parser.add_argument("--systems-top-n", type=int, default=15)
    parser.add_argument("--systems-max-tf", type=int, default=46)
    parser.add_argument("--systems-fdr", type=float, default=0.05)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--confirm-heavy", action="store_true")
    parser.add_argument("--plan-only", action="store_true")
    parser.add_argument("--reset-paths", action="store_true")
    return parser


class Runner:
    def __init__(self, args: argparse.Namespace):
        self.args = args
        self.profile_file = profile_path(args.path_config)
        if args.reset_paths and self.profile_file.exists():
            self.profile_file.unlink()
        self.profile = load_profile(self.profile_file)
        self.paths = self.resolve_paths()
        self.analysis = Path(self.paths["analysis_root"]) / args.analysis_project
        self.logs = self.analysis / "00_logs"
        self.prs = self.analysis / "13_prs"
        self.rscript = ""
        self.model_python = ""

    def resolve_paths(self) -> dict[str, str]:
        ## Huang-lab public layout is the default. The legacy UKB_data tree is
        ## used only when the user explicitly supplies/saves it as DIR0.
        base = local_path(self.args.dir0 or self.profile.get("dir0", "") or "D:/")
        legacy = Path(base).name.lower() == "ukb_data"
        if legacy and Path(base).name.lower() == "ukb_data":
            defaults = {
                "phe_dir": join_root(base, "phe"), "analysis_root": join_root(base, "analysis"),
                "raw_protein_file": join_root(base, "phe/raw/prot_full_unimputed.tsv"),
                "phenotype_rds": join_root(base, "phe/Rdata/all.rds"),
                "raw_phenotype_file": join_root(base, "pheno.tsv.gz"),
                "panel_mapping_file": join_root(base, "ppp/map.raw/olink_protein_map_3k_v1.tsv"),
                "cmr_feature_file": join_root(base, "analysis/sleepchart_reproduction/data/mribag_features/heart/feature.tsv"),
                "pqtl_root": join_root(base, "ppp/clean"),
            }
        else:
            defaults = {
                "phe_dir": join_root(base, "data/ukb/phe"), "analysis_root": join_root(base, "analysis"),
                "raw_protein_file": join_root(base, "data/ukb/phe/raw/prot_full_unimputed.tsv"),
                "phenotype_rds": join_root(base, "data/ukb/phe/Rdata/all.rds"),
                "raw_phenotype_file": join_root(base, "data/ukb/phe/pheno.tsv.gz"),
                "panel_mapping_file": join_root(base, "data.BIG/gwas/ppp/olink_protein_map_3k_v1.tsv"),
                "cmr_feature_file": join_root(base, "analysis/sleepchart_reproduction/data/mribag_features/heart/feature.tsv"),
                "pqtl_root": join_root(base, "ppp/clean"),
            }
        cli = {
            "phe_dir": self.args.phe_dir, "analysis_root": self.args.analysis_root,
            "raw_protein_file": self.args.raw_protein_file, "phenotype_rds": self.args.phenotype_rds,
            "raw_phenotype_file": self.args.raw_phenotype_file,
            "panel_mapping_file": self.args.panel_mapping_file,
            "cmr_feature_file": self.args.cmr_feature_file, "pqtl_root": self.args.pqtl_root,
        }
        result = {"dir0": base}
        for key, default in defaults.items():
            result[key] = local_path(cli[key] or self.profile.get(key, "") or default)
        result.update({
            "supplement_workbook_file": local_path(self.args.supplement_workbook_file or self.profile.get("supplement_workbook_file", "") or str(PROJECT / "references/raw/pwaf072_supplementary_table_1.xlsx")),
            "supplement_methods_file": local_path(self.args.supplement_methods_file or self.profile.get("supplement_methods_file", "") or str(PROJECT / "references/raw/pwaf072_supplementary_figure_1.pdf")),
            "olink_processing_start_date_file": local_path(self.args.olink_processing_start_date_file or self.profile.get("olink_processing_start_date_file", "") or str(PROJECT / "references/raw/olink_processing_start_date.dat")),
            "genotype_root": self.args.genotype_root or self.profile.get("genotype_root", "") or "Z:/projects/genotype_pc_nas/imputed_pgen_autosomes",
        })
        result["mr_outcome_lookup_dir"] = local_path(
            self.args.mr_outcome_lookup_dir or self.profile.get("mr_outcome_lookup_dir", "")
            or join_root(result["analysis_root"], self.args.analysis_project, "11_mr/outcome_lookup")
        )
        return result

    def save_paths(self) -> None:
        values = {key: value for key, value in self.paths.items() if value}
        save_profile(self.profile_file, values)
        print(f"Saved paths: {self.profile_file}")

    def require_path(self, key: str, label: str, kind: str = "file") -> None:
        value = self.paths[key]
        valid = existing_file(value) if kind == "file" else existing_dir(value)
        if valid:
            return
        if self.args.plan_only:
            print(f"WARN missing {label}: {value}")
            return
        if self.args.path_prompt == "off" or not sys.stdin.isatty():
            fail(f"{label} not found: {value}. Run './yu.sh setup' or provide its option.")
        while True:
            answer = input(f"{label}\nCurrent: {value}\nEnter path (Q to cancel): ").strip()
            if answer.lower() in ("q", "quit", "exit"):
                fail(f"Path setup cancelled: {label}")
            answer = local_path(answer)
            valid = existing_file(answer) if kind == "file" else existing_dir(answer)
            if valid:
                self.paths[key] = answer
                return
            print(f"Not found: {answer}")

    def resolve_inputs(self, steps: list[int]) -> None:
        if any(step in steps for step in (1, 2, 3, 4, 5, 7, 9)):
            self.require_path("raw_protein_file", "Unimputed protein table")
            self.require_path("phenotype_rds", "Phenotype RDS")
        if any(step in steps for step in (1, 2, 3, 4, 5, 7)):
            self.require_path("raw_phenotype_file", "Raw phenotype file")
            self.require_path("panel_mapping_file", "Olink panel mapping")
        if any(step in steps for step in (1, 7, 8, 11)):
            self.require_path("supplement_workbook_file", "Official supplementary workbook")
        if 1 in steps:
            self.require_path("supplement_methods_file", "Official supplementary methods")
        if any(step in steps for step in (1, 2)):
            self.require_path("olink_processing_start_date_file", "Olink processing-date table")
        if 5 in steps:
            self.require_path("cmr_feature_file", "CMR feature table")
        if 6 in steps:
            self.require_path("pqtl_root", "pQTL root", "dir")
        if self.args.protein_panel == "custom" and not self.args.model_proteins:
            if not self.args.model_protein_file:
                fail("--protein-panel custom requires --model-protein-file or --model-proteins.")
            if not existing_file(self.args.model_protein_file):
                fail(f"Custom protein file not found: {self.args.model_protein_file}")

    def show_plan(self, steps: list[int]) -> None:
        print("\nResolved paths")
        for key, value in self.paths.items():
            print(f"  {key:34} {value}")
        print(f"  {'project':34} {PROJECT}")
        print(f"  {'analysis':34} {self.analysis}")
        print("\nExecution")
        for step in steps:
            print(f"  {step:>2}  {STEPS[step][0]:42} {STEPS[step][1]}")
        print(f"  disease={self.args.disease}; protein_panel={self.args.protein_panel}")

    def ensure_r(self) -> None:
        if not self.rscript:
            self.rscript = resolve_rscript(self.args.rscript)

    def ensure_python(self) -> None:
        if not self.model_python:
            self.model_python = resolve_model_python(self.args.model_python)
        code = "import sys,lightgbm;print('.'.join(map(str,sys.version_info[:3])),lightgbm.__version__)"
        result = subprocess.run([self.model_python, "-c", code], text=True, capture_output=True)
        if result.returncode or not valid_model_runtime(result.stdout):
            fail(f"Formal models require Python 3.9 and LightGBM 3.3.2. Found: {result.stdout.strip() or result.stderr.strip()}")

    def r_arguments(self, stage: str, overrides: dict[str, object] | None = None, entry: Path = FULL_R) -> list[str]:
        self.ensure_r()
        values: dict[str, object] = {
            "mode": stage, "dir0": self.paths["dir0"], "analysis_root": self.paths["analysis_root"],
            "analysis_project": self.args.analysis_project, "raw_protein_file": self.paths["raw_protein_file"],
            "cmr_feature_file": self.paths["cmr_feature_file"], "cmr_metric_subset": "all",
            "phenotype_rds": self.paths["phenotype_rds"], "raw_phenotype_file": self.paths["raw_phenotype_file"],
            "panel_mapping_file": self.paths["panel_mapping_file"],
            "supplement_workbook_file": self.paths["supplement_workbook_file"],
            "supplement_methods_file": self.paths["supplement_methods_file"],
            "olink_processing_start_date_file": self.paths["olink_processing_start_date_file"],
            "endpoint_subset": self.args.disease, "workers": self.args.workers,
            "bootstrap_n": self.args.bootstrap_n,
            "prediction_panel_mode": self.args.protein_panel,
            "figure4_extra_projects": self.args.figure4_extra_project,
            "figure4_extra_outcomes": self.args.figure4_extra_outcome,
            "figure4_extra_labels": self.args.figure4_extra_label,
            "pqtl_root": self.paths["pqtl_root"], "mr_outcome_lookup_dir": self.paths["mr_outcome_lookup_dir"],
            "cmest_shard_index": 1, "cmest_shard_count": 1,
            "cmest_pilot_nboot": self.args.cmest_pilot_boot,
            "string_required_score": self.args.systems_score,
            "systems_top_n_per_outcome": self.args.systems_top_n,
            "systems_max_tf": self.args.systems_max_tf,
            "systems_enrichment_fdr": self.args.systems_fdr,
        }
        if entry == PRS_R:
            values = {key: values[key] for key in ("mode", "dir0", "analysis_root", "analysis_project", "raw_protein_file", "phenotype_rds", "endpoint_subset", "workers")}
        if overrides:
            values.update(overrides)
        path_keys = {
            key for key in values
            if key == "dir0" or key.endswith("_file") or key.endswith("_root") or key.endswith("_dir")
        }
        arguments = [self.rscript, "--vanilla", runtime_path(entry, self.rscript)]
        for key, value in values.items():
            if key in path_keys and value:
                value = runtime_path(str(value), self.rscript)
            arguments.append(f"--{key}={value}")
        if self.args.resume:
            arguments.append("--resume=true")
        if self.args.force:
            arguments.append("--force=true")
        return arguments

    def run_r(self, stage: str, overrides: dict[str, object] | None = None, entry: Path = FULL_R) -> None:
        self.logs.mkdir(parents=True, exist_ok=True)
        log = self.logs / f"{stage}_{datetime.now():%Y%m%d_%H%M%S}.log"
        print(f"{stamp()} | R stage={stage}")
        tee_run(self.r_arguments(stage, overrides, entry), log)

    def python_arguments(self, stage: str) -> list[str]:
        self.ensure_python()
        values = [
            self.model_python, runtime_path(MODEL_PY, self.model_python), "--mode", stage,
            "--project-dir", runtime_path(PROJECT, self.model_python),
            "--analysis-dir", runtime_path(self.analysis, self.model_python),
            "--raw-protein-file", runtime_path(self.paths["raw_protein_file"], self.model_python),
            "--panel-mapping-file", runtime_path(self.paths["panel_mapping_file"], self.model_python),
            "--config", runtime_path(MODEL_CONFIG, self.model_python),
            "--endpoint-subset", self.args.disease, "--workers", str(self.args.workers),
            "--model-jobs", str(self.args.model_jobs),
            "--prediction-panel-mode", self.args.protein_panel,
        ]
        if self.args.model_protein_file:
            values.extend(["--custom-protein-panel-file", runtime_path(local_path(self.args.model_protein_file), self.model_python)])
        if self.args.model_proteins:
            values.extend(["--custom-proteins", self.args.model_proteins])
        if self.args.resume:
            values.append("--resume")
        return values

    def run_python(self, stage: str) -> None:
        log = self.logs / f"{stage}_{datetime.now():%Y%m%d_%H%M%S}.log"
        print(f"{stamp()} | Python stage={stage}")
        tee_run(self.python_arguments(stage), log)

    def requested_endpoints(self) -> list[str]:
        all_ids = read_first_column(F_DIR / "config" / "outcomes.csv")
        if self.args.disease.lower() == "all":
            return all_ids
        requested = list(dict.fromkeys(item.strip() for item in self.args.disease.split(",") if item.strip()))
        unknown = sorted(set(requested) - set(all_ids))
        if unknown:
            fail(f"Unknown disease IDs: {', '.join(unknown)}")
        return requested

    def parallel_r(
        self, label: str, items: list[str], jobs: int, stage: str, overrides, validate,
        workers_per_job: int | None = None,
    ) -> None:
        jobs = min(max(1, jobs), len(items))
        shard_workers = workers_per_job or max(1, self.args.workers // jobs)
        print(f"{stamp()} | {label} items={len(items)} concurrent={jobs} workers/job={shard_workers}")

        def task(item: str) -> tuple[str, str]:
            stdout = self.logs / f"{label.lower()}_{item}_{datetime.now():%Y%m%d_%H%M%S}.out.log"
            stderr = stdout.with_name(stdout.name.replace(".out.log", ".err.log"))
            command = self.r_arguments(stage, {**overrides(item), "workers": shard_workers})
            started = time.monotonic()
            with stdout.open("w", encoding="utf-8") as out, stderr.open("w", encoding="utf-8") as err:
                result = subprocess.run(command, stdout=out, stderr=err)
            if result.returncode or not validate(item):
                return item, f"exit={result.returncode}; stderr={stderr}"
            return item, f"PASS minutes={(time.monotonic() - started) / 60:.1f}"

        failures = []
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            futures = {pool.submit(task, item): item for item in items}
            for future in as_completed(futures):
                item, status = future.result()
                print(f"{stamp()} | {label} {item} | {status}", flush=True)
                if not status.startswith("PASS"):
                    failures.append(item)
        if failures:
            fail(f"{label} shards failed: {', '.join(failures)}. Completed shards remain reusable with --resume.")

    def cox(self) -> None:
        self.run_r("cox_prepare")
        def valid(endpoint: str) -> bool:
            base = self.analysis / "05_cox"
            expected = [
                base / f"full_incident_{endpoint}_cox.csv.gz", base / f"full_incident_{endpoint}_cox.contract.json",
                base / f"derivation_{endpoint}_cox.csv.gz", base / f"derivation_{endpoint}_cox.contract.json",
            ]
            marker = self.logs / f"cox_shard_{re.sub('[^A-Za-z0-9]+', '_', endpoint)}.done.json"
            try:
                marker_ok = json.loads(marker.read_text(encoding="utf-8-sig")).get("status") == "PASS"
            except Exception:
                marker_ok = False
            return marker_ok and all(path.is_file() and path.stat().st_size for path in expected)
        self.parallel_r("Cox", self.requested_endpoints(), self.args.cox_jobs, "cox_shard", lambda x: {"endpoint_subset": x}, valid)
        self.run_r("cox_merge")

    def cmr(self) -> None:
        self.run_r("cmr_prepare")
        metrics = read_first_column(F_DIR / "config" / "cmr_metrics.csv")
        base = self.analysis / "07_cmr"
        valid = lambda item: all((base / f"cmr_{item}.{suffix}").is_file() for suffix in ("csv.gz", "contract.json"))
        self.parallel_r("CMR", metrics, self.args.cmr_jobs, "cmr_shard", lambda x: {"cmr_metric_subset": x}, valid)
        self.run_r("cmr_merge")

    def cmest(self) -> None:
        self.run_r("mediation_cmest_pilot")
        count = self.args.cmest_jobs
        base = self.analysis / "12_mediation" / "cmest_shards"
        items = [str(index) for index in range(1, count + 1)]
        def valid(item: str) -> bool:
            prefix = base / f"cmest_shard_{int(item):03d}_of_{count:03d}"
            return Path(f"{prefix}.csv").is_file() and Path(f"{prefix}.contract.json").is_file()
        self.parallel_r(
            "CMAverse", items, count, "mediation_cmest_shard",
            lambda x: {"cmest_shard_index": int(x), "cmest_shard_count": count}, valid,
            workers_per_job=1,
        )
        self.run_r("mediation_cmest_merge", {"cmest_shard_count": count})

    def prepare_gwas(self) -> None:
        command = [
            sys.executable, str(GWAS_PREP), "--project-dir", str(F_DIR),
            "--analysis-dir", str(self.analysis), "--max-p", "0.0005",
        ]
        if self.args.resume:
            command.append("--resume")
        tee_run(command, self.logs / f"prs_prepare_gwas_{datetime.now():%Y%m%d_%H%M%S}.log")

    def install_plink2(self) -> str:
        self.ensure_python()
        windows = self.model_python.lower().endswith(".exe")
        if not windows:
            found = shutil.which("plink2")
            if found:
                return found
            fail("PLINK2 not found for the active Python runtime.")
        target = F_DIR / "tools" / "bin" / "plink2.exe"
        if target.is_file():
            return str(target)
        target.parent.mkdir(parents=True, exist_ok=True)
        url = "https://s3.amazonaws.com/plink2-assets/alpha6/plink2_win_avx2_20250129.zip"
        archive = target.parent / "plink2.zip"
        print(f"{stamp()} | Installing PLINK2 (~8 MB)")
        urllib.request.urlretrieve(url, archive)
        with zipfile.ZipFile(archive) as package:
            member = next(name for name in package.namelist() if name.lower().endswith("plink2.exe"))
            with package.open(member) as source, target.open("wb") as destination:
                shutil.copyfileobj(source, destination)
        archive.unlink()
        return str(target)

    def score_prs(self) -> None:
        self.ensure_python()
        plink2 = self.install_plink2()
        command = [
            self.model_python, runtime_path(PRS_SCORER, self.model_python),
            "--project-dir", runtime_path(F_DIR, self.model_python),
            "--analysis-dir", runtime_path(self.analysis, self.model_python),
            "--nas-root", runtime_path(self.paths["genotype_root"], self.model_python),
            "--plink2", runtime_path(plink2, self.model_python),
            "--workers", str(self.args.workers), "--score-jobs", str(self.args.score_jobs),
            "--memory-mb", str(self.args.memory_mb), "--start-chr", "1", "--end-chr", "22",
        ]
        if self.args.resume:
            command.append("--resume")
        tee_run(command, self.logs / f"prs_score_{datetime.now():%Y%m%d_%H%M%S}.log")

    def prs_associations(self) -> None:
        endpoints = read_first_column(F_DIR / "config" / "prs_gwas_sources.tsv", "\t")
        base = self.prs
        def valid(endpoint: str) -> bool:
            result = base / f"prs_protein_associations_{endpoint}.csv.gz"
            marker = base / f"prs_association_{endpoint}.done.json"
            try:
                return result.is_file() and result.stat().st_size > 0 and json.loads(marker.read_text(encoding="utf-8-sig")).get("status") == "PASS"
            except Exception:
                return False
        self.parallel_r("PRS", endpoints, self.args.association_jobs, "associate_shard", lambda x: {"endpoint_subset": x}, valid)

    def run_prs(self) -> None:
        self.run_r("preflight", {"endpoint_subset": "all"}, PRS_R)
        self.prepare_gwas()
        self.score_prs()
        self.run_r("merge_scores", {"endpoint_subset": "all"}, PRS_R)
        self.prs_associations()
        self.run_r("merge_associations", {"endpoint_subset": "all"}, PRS_R)
        self.run_r("figures", {"endpoint_subset": "all", "prediction_panel_mode": "local_reselected"})
        self.run_r("report", {"endpoint_subset": "all"}, PRS_R)

    def systems(self) -> None:
        self.ensure_r()
        installer = F_DIR / "tools" / "install_yu_systems_packages.R"
        tee_run(
            [self.rscript, "--vanilla", runtime_path(installer, self.rscript)],
            self.logs / f"systems_packages_{datetime.now():%Y%m%d_%H%M%S}.log",
        )
        for stage in ("systems_prepare", "systems_enrichment", "systems_tf", "systems_ppi", "systems_figures"):
            self.run_r(stage)

    def assert_final_inputs(self) -> None:
        required = [
            self.analysis / "14_enrichment/enrichment_results.csv",
            self.analysis / "14_enrichment/local_cox_systems/figure6c_cvd_protein_tf_edges.csv",
            self.analysis / "14_enrichment/local_cox_systems/figure6d_main_cluster_nodes.csv",
            self.analysis / "14_enrichment/local_cox_systems/figure6d_main_cluster_edges.csv.gz",
            self.analysis / "15_figures/figure6a_local_prs_protein_heatmap.png",
        ]
        missing = [str(path) for path in required if not path.is_file() or not path.stat().st_size]
        if missing:
            fail("Final figure inputs are incomplete. Run Step 10 after Steps 1-9. Missing:\n" + "\n".join(missing))

    def assert_deliverables(self) -> None:
        prefixes = [
            "figure1_workflow", "figure2_incident_protein_associations", "figure3_local_cmr",
            "figure4_prediction_and_importance", "figure5_local_mr_and_mediation",
            "figure6abcd_local_systems_biology",
        ]
        required = [self.analysis / "15_figures" / f"{prefix}.{ext}" for prefix in prefixes for ext in ("pdf", "png", "tiff")]
        required.append(self.analysis / "17_report/RESULTS_AND_QC.md")
        missing = [str(path) for path in required if not path.is_file() or not path.stat().st_size]
        if missing:
            fail("Deliverable validation failed:\n" + "\n".join(missing))
        print("FINAL PACKAGE PASS | Figures 1-6: PDF, PNG, TIFF; report: PASS")

    def step(self, number: int) -> None:
        if number == 1:
            self.run_r("sources"); self.run_r("preflight"); self.ensure_python()
        elif number == 2:
            self.run_r("cohort")
        elif number == 3:
            self.cox()
        elif number == 4:
            self.run_python("select"); self.run_python("train"); self.run_r("evaluate")
        elif number == 5:
            self.cmr()
        elif number == 6:
            self.run_r("mr_prepare"); self.run_r("mr_run")
        elif number == 7:
            self.run_r("mediation_prepare"); self.run_r("mediation_run")
        elif number == 8:
            self.cmest()
        elif number == 9:
            self.run_prs()
        elif number == 10:
            self.systems()
        elif number == 11:
            self.assert_final_inputs(); self.run_r("figures"); self.run_r("systems_figures"); self.run_r("report"); self.assert_deliverables()

    def doctor(self) -> None:
        prediction_only = os.environ.get("YU_DOCTOR_SCOPE", "").strip().lower() == "prediction"
        required = {
            "R entry": FULL_R, "Python model": MODEL_PY, "protein": Path(self.paths["raw_protein_file"]),
            "phenotype": Path(self.paths["phenotype_rds"]), "raw phenotype": Path(self.paths["raw_phenotype_file"]),
            "panel mapping": Path(self.paths["panel_mapping_file"]),
            "Olink dates": Path(self.paths["olink_processing_start_date_file"]),
        }
        if not prediction_only:
            required["supplement tables"] = Path(self.paths["supplement_workbook_file"])
            required["supplement methods"] = Path(self.paths["supplement_methods_file"])
        missing = []
        for label, path in required.items():
            ok = path.is_file() and path.stat().st_size > 0
            print(f"{label:24} {'PASS' if ok else 'FAIL'}  {path}")
            if not ok:
                missing.append(str(path))
        self.ensure_r(); self.ensure_python()
        package_specification = (
            "data.table,R.utils,jsonlite,digest,readxl,survival,pROC,ggplot2,bit64"
            if prediction_only else
            "data.table,R.utils,jsonlite,digest,readxl,survival,pROC,ggplot2,bit64,patchwork,ggrepel,ragg,httr2,curl,igraph,ggalluvial,scales,msigdbr,AnnotationDbi,org.Hs.eg.db,GO.db,TwoSampleMR,CMAverse"
        )
        packages = package_specification.split(",")
        expression = "p=c(" + ",".join(repr(item) for item in packages) + ");m=p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)];cat(if(length(m))paste(m,collapse=',')else'PASS')"
        result = subprocess.run([self.rscript, "--vanilla", "-e", expression], text=True, capture_output=True)
        status = result.stdout.strip().splitlines()[-1] if result.stdout.strip() else result.stderr.strip()
        package_label = "R packages (steps 1-4)" if prediction_only else "R packages"
        print(f"{package_label:24} {'PASS' if result.returncode == 0 and status == 'PASS' else 'FAIL'}  {status}")
        if result.returncode or status != "PASS":
            missing.append(f"R packages: {status}")
        if not prediction_only:
            if required["supplement tables"].is_file() and sha256(required["supplement tables"]) != EXPECTED_XLSX:
                missing.append("supplement workbook SHA256")
            if required["supplement methods"].is_file() and sha256(required["supplement methods"]) != EXPECTED_PDF:
                missing.append("supplement methods SHA256")
        if missing:
            fail("Doctor failed:\n" + "\n".join(missing))
        print("YU PROJECT DOCTOR PASS")

    def status(self) -> None:
        print(f"Analysis: {self.analysis}")
        try:
            processes = subprocess.check_output(["ps", "-eo", "pid,lstart,args"], text=True)
            active = [line for line in processes.splitlines() if any(token in line for token in ("99_run_yu", "04_full_reproduction.py", "score_prs_directnas", "plink2")) and "yu_runner.py status" not in line]
            print("\nActive processes")
            print("\n".join(active) if active else "none")
        except Exception:
            pass
        print("\nStage markers")
        markers = sorted(self.logs.glob("*.done.json"), key=lambda path: path.stat().st_mtime) if self.logs.is_dir() else []
        for marker in markers[-30:]:
            print(f"{datetime.fromtimestamp(marker.stat().st_mtime):%Y-%m-%d %H:%M:%S}  {marker.name}")
        run_log = self.logs / "run.log"
        if run_log.is_file():
            print("\nLatest run.log")
            print("".join(run_log.read_text(encoding="utf-8", errors="replace").splitlines(True)[-40:]))
        for name in ("gwas_prepare_status.tsv", "score_status.tsv", "score_variant_counts.tsv", "prs_association_summary.csv"):
            path = self.prs / name
            if path.is_file():
                print(f"\n{name}")
                print("".join(path.read_text(encoding="utf-8", errors="replace").splitlines(True)[-12:]))

    def install(self) -> None:
        self.ensure_r()
        try:
            python = resolve_model_python(self.args.model_python)
        except RuntimeError:
            conda = resolve_conda(self.args.conda)
            subprocess.run([conda, "create", "-y", "-n", "yu_proteomic_repo_py39", "python=3.9"], check=True)
            result = subprocess.check_output([conda, "env", "list", "--json"], text=True)
            env_dir = next(path for path in json.loads(result)["envs"] if path.endswith("yu_proteomic_repo_py39"))
            python = str(Path(env_dir) / ("python.exe" if os.name == "nt" or conda.lower().endswith(".exe") else "bin/python"))
        requirements = runtime_path(F_DIR / "config/requirements-py39.txt", python)
        subprocess.run([python, "-m", "pip", "install", "-r", requirements], check=True)
        manifest = PROJECT / "dist" / "dependency_manifests"
        manifest.mkdir(parents=True, exist_ok=True)
        subprocess.run([
            self.rscript, "--vanilla", runtime_path(F_DIR / "tools/install_yu_r_dependencies.R", self.rscript),
            runtime_path(PROJECT, self.rscript), runtime_path(manifest / "r_dependencies.txt", self.rscript),
        ], check=True)
        print(f"INSTALL PASS\nRscript: {self.rscript}\nPython: {python}")

    def package(self) -> None:
        output = PROJECT / "dist"
        output.mkdir(exist_ok=True)
        name = f"Yu_protein_analysis_code_{datetime.now():%Y%m%d}.zip"
        target = output / name
        tracked = subprocess.check_output(
            ["git", "ls-files", "--cached", "--others", "--exclude-standard"],
            cwd=PROJECT, text=True,
        ).splitlines()
        with tempfile.TemporaryDirectory() as temp:
            stage = Path(temp) / "yu-protein-analysis"
            for relative in tracked:
                source = PROJECT / relative
                if not source.is_file():
                    continue
                destination = stage / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, destination)
            manifest_rows = []
            for path in sorted(item for item in stage.rglob("*") if item.is_file()):
                manifest_rows.append((path.relative_to(stage).as_posix(), path.stat().st_size, sha256(path)))
            manifest = stage / "FILE_MANIFEST_SHA256.csv"
            with manifest.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["relative_path", "size_bytes", "sha256"])
                writer.writerows(manifest_rows)
            if target.exists():
                target.unlink()
            with zipfile.ZipFile(target, "w", zipfile.ZIP_DEFLATED) as archive:
                for path in stage.rglob("*"):
                    if path.is_file():
                        archive.write(path, path.relative_to(stage.parent))
        digest = sha256(target)
        target.with_suffix(target.suffix + ".sha256").write_text(f"{digest}  {target.name}\n", encoding="ascii")
        print(f"PACKAGE PASS\nZIP: {target}\nSHA256: {digest}")


def main() -> None:
    args = build_parser().parse_args()
    command = (args.step or args.command).lower()
    if command in ("help", "-h", "--help"):
        build_parser().print_help()
        return
    runner = Runner(args)
    if command == "paths":
        runner.show_plan([])
        print(f"\nProfile: {runner.profile_file}")
        return
    if command == "setup":
        runner.resolve_inputs([1, 2, 3, 4])
        runner.save_paths()
        runner.show_plan([])
        print("Path setup complete. No analysis started.")
        return
    if command == "status":
        runner.status(); return
    if command == "package":
        runner.package(); return
    if command == "install":
        runner.install(); return
    steps = [1, 2, 3, 4] if command == "doctor" else parse_steps(command)
    runner.resolve_inputs(steps)
    runner.show_plan(steps)
    if args.plan_only:
        return
    runner.save_paths()
    if command == "doctor":
        runner.doctor(); return
    if (args.disease != "all" or args.protein_panel != "local_reselected") and args.analysis_project == "yu_proteomic_repo_v3":
        fail("yu_proteomic_repo_v3 is protected. Use a new --analysis-project for alternate disease or protein panels.")
    if any(step in (8, 9) for step in steps) and not args.confirm_heavy:
        fail("Steps 8 and 9 are long-running. Review with --plan-only, then add --confirm-heavy.")
    runner.logs.mkdir(parents=True, exist_ok=True)
    state_file = runner.logs / f"step_runner_{datetime.now():%Y%m%d_%H%M%S}.json"
    state = {"status": "RUNNING", "steps": steps, "completed": [], "started_at": stamp(), "analysis": str(runner.analysis)}
    state_file.write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")
    try:
        for step in steps:
            print(f"{stamp()} | STEP {step} START | {STEPS[step][0]}")
            runner.step(step)
            state["completed"].append(step)
            state_file.write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")
            print(f"{stamp()} | STEP {step} DONE")
        state.update(status="PASS", ended_at=stamp())
    except Exception as error:
        state.update(status="ERROR", ended_at=stamp(), error=str(error))
        raise
    finally:
        state_file.write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")
    print(f"{stamp()} | COMPLETE | steps={','.join(map(str, steps))} | analysis={runner.analysis}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr)
        raise SystemExit(130)
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
