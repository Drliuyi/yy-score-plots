#!/usr/bin/env python3
"""Fit one locked outer fold of the all-protein Yu fair comparison."""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from pathlib import Path

import lightgbm as lgb
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=["fold", "status"], required=True)
    parser.add_argument("--fold", type=int)
    parser.add_argument("--common-root", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--workers", type=int, default=8)
    return parser.parse_args()


def atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = ".csv.gz" if str(path).endswith(".gz") else ".csv"
    temporary = path.with_name(f"{path.stem}.tmp.{os.getpid()}{suffix}")
    frame.to_csv(temporary, index=False)
    os.replace(temporary, path)


def atomic_text(lines: list[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    temporary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def load_matrix(path: Path, rows: int, columns: int) -> np.memmap:
    expected = rows * columns * 4
    if not path.exists() or path.stat().st_size != expected:
        raise RuntimeError(f"Float32 matrix contract failed: {path}")
    raw = np.memmap(path, dtype="<f4", mode="r", shape=(rows * columns,))
    return raw.reshape((rows, columns), order="F")


def known_horizon(frame: pd.DataFrame, horizon: float) -> tuple[np.ndarray, np.ndarray]:
    time_value = frame["time"].to_numpy(dtype=float)
    event = frame["event"].to_numpy(dtype=int)
    known = ((event == 1) & (time_value <= horizon)) | (time_value > horizon)
    label = ((event == 1) & (time_value <= horizon)).astype(np.int8)
    return known, label


def logit(probability: np.ndarray) -> np.ndarray:
    value = np.clip(probability, 1e-7, 1 - 1e-7)
    return np.log(value / (1 - value))


def main() -> None:
    args = parse_args()
    common = Path(args.common_root)
    output = Path(args.output_root)
    with open(args.config, "r", encoding="utf-8") as handle:
        cfg = json.load(handle)
    expected = cfg["expected"]
    markers = [output / "models" / f"fold{fold:02d}_COMPLETE" for fold in range(1, 6)]
    if args.stage == "status":
        for fold, marker in enumerate(markers, start=1):
            print(f"fold{fold:02d}\t{'COMPLETE' if marker.exists() else 'MISSING'}")
        raise SystemExit(0 if all(path.exists() for path in markers) else 1)
    if args.fold is None or args.fold not in range(1, 6):
        raise ValueError("--fold=1..5 is required")
    fold = args.fold
    required = [
        common / "COMPLETE",
        common / "participants_yin.csv.gz",
        common / "participants_yang.csv.gz",
        common / "protein_features.csv",
        common / "protein_yin.f32",
        common / "protein_yang.f32",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise RuntimeError("Fair common input missing: " + "; ".join(missing))
    prediction_path = output / "predictions" / f"fold{fold:02d}_yin.csv.gz"
    yang_path = output / "predictions" / f"fold{fold:02d}_yang.csv.gz"
    model_path = output / "models" / f"fold{fold:02d}_model.txt"
    diagnostic_path = output / "models" / f"fold{fold:02d}_diagnostic.json"
    marker = output / "models" / f"fold{fold:02d}_COMPLETE"
    fold_outputs = [prediction_path, yang_path, model_path, diagnostic_path, marker]
    if all(path.exists() for path in fold_outputs):
        print(f"YU_FAIR_FOLD_PRESERVED fold={fold}")
        return
    if any(path.exists() for path in fold_outputs):
        raise RuntimeError(f"Partial Yu fair fold exists; refusing overwrite: fold {fold}")

    yin = pd.read_csv(common / "participants_yin.csv.gz", dtype={"eid": str})
    yang = pd.read_csv(common / "participants_yang.csv.gz", dtype={"eid": str})
    features = pd.read_csv(common / "protein_features.csv")["feature"].astype(str).tolist()
    if (
        len(yin) != int(expected["yin_n"])
        or int(yin["event"].sum()) != int(expected["yin_events"])
        or len(yang) != int(expected["yang_n"])
        or len(features) != int(expected["protein_n"])
        or yin["eid"].duplicated().any()
        or yang["eid"].duplicated().any()
    ):
        raise RuntimeError("Yu fair participant/feature contract failed")
    x_yin = load_matrix(common / "protein_yin.f32", len(yin), len(features))
    x_yang = load_matrix(common / "protein_yang.f32", len(yang), len(features))
    train_outer = yin["outer_fold"].to_numpy(dtype=int) != fold
    test_index = np.flatnonzero(~train_outer)
    train_index_all = np.flatnonzero(train_outer)
    known, label_all = known_horizon(yin.iloc[train_index_all], float(cfg["horizon_years"]))
    train_index = train_index_all[known]
    y_train = label_all[known]
    if np.unique(y_train).size != 2:
        raise RuntimeError("Yu fair training label has fewer than two classes")

    parameters = cfg["yu_fair"]
    model = lgb.LGBMClassifier(
        objective=str(parameters["objective"]),
        n_estimators=int(parameters["n_estimators"]),
        max_depth=int(parameters["max_depth"]),
        num_leaves=int(parameters["num_leaves"]),
        subsample=float(parameters["subsample"]),
        subsample_freq=int(parameters["subsample_freq"]),
        learning_rate=float(parameters["learning_rate"]),
        colsample_bytree=float(parameters["colsample_bytree"]),
        max_bin=int(parameters["max_bin"]),
        device_type=str(parameters["device_type"]),
        random_state=int(cfg["seed"]) + fold * 1000,
        n_jobs=max(1, int(args.workers)),
        verbosity=-1,
    )
    print(
        f"YU_FAIR_FIT fold={fold} train_known_n={len(train_index)} "
        f"cases={int(y_train.sum())} test_n={len(test_index)} device={parameters['device_type']}"
    )
    started = time.time()
    model.fit(np.asarray(x_yin[train_index, :]), y_train, feature_name=features)
    yin_probability = model.predict_proba(np.asarray(x_yin[test_index, :]))[:, 1]
    yang_probability = model.predict_proba(np.asarray(x_yang))[:, 1]
    yin_score = logit(yin_probability)
    yang_score = logit(yang_probability)
    if not np.isfinite(yin_score).all() or not np.isfinite(yang_score).all():
        raise RuntimeError("Yu fair produced non-finite scores")

    yin_output = yin.iloc[test_index][["eid", "outer_fold", "time", "event"]].copy()
    yin_output["method"] = "yu-fair"
    yin_output["score_raw"] = yin_score
    yang_output = pd.DataFrame(
        {
            "eid": yang["eid"].astype(str),
            "outer_fold": fold,
            "method": "yu-fair",
            "score_raw": yang_score,
        }
    )
    model_path.parent.mkdir(parents=True, exist_ok=True)
    stage_model = model_path.with_name(f"{model_path.name}.tmp.{os.getpid()}")
    model.booster_.save_model(str(stage_model))
    diagnostic = {
        "method": "yu-fair",
        "outer_fold": fold,
        "target": "five-year incident CAD",
        "training_n": int(len(train_index)),
        "training_cases": int(y_train.sum()),
        "test_n": int(len(test_index)),
        "feature_n": len(features),
        "elapsed_seconds": time.time() - started,
        "lightgbm_version": lgb.__version__,
        "parameters": parameters,
    }
    stage_diagnostic = diagnostic_path.with_name(f"{diagnostic_path.name}.tmp.{os.getpid()}")
    stage_diagnostic.write_text(json.dumps(diagnostic, indent=2) + "\n", encoding="utf-8")
    atomic_csv(yin_output, prediction_path)
    atomic_csv(yang_output, yang_path)
    os.replace(stage_model, model_path)
    os.replace(stage_diagnostic, diagnostic_path)
    atomic_text(
        [
            "status=COMPLETE",
            f"fold={fold}",
            f"training_n={len(train_index)}",
            f"training_cases={int(y_train.sum())}",
            f"completed_at={time.strftime('%Y-%m-%dT%H:%M:%S%z')}",
        ],
        marker,
    )
    print(f"YU_FAIR_FOLD_COMPLETE fold={fold}")


if __name__ == "__main__":
    main()
