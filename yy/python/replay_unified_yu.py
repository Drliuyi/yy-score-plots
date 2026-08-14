#!/usr/bin/env python3
"""Replay five frozen Yu-style outer models without fitting or tuning."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import lightgbm as lgb
import numpy as np
import pandas as pd


def logit(probability: np.ndarray) -> np.ndarray:
    value = np.clip(np.asarray(probability, dtype=float), 1e-7, 1 - 1e-7)
    return np.log(value / (1 - value))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-root", required=True)
    parser.add_argument("--benchmark-root", required=True)
    parser.add_argument("--yin-output", required=True)
    parser.add_argument("--yang-output", required=True)
    parser.add_argument("--metadata-output", required=True)
    parser.add_argument("--qc-output", required=True)
    args = parser.parse_args()

    cache = Path(args.cache_root)
    benchmark = Path(args.benchmark_root)
    yin = pd.read_csv(cache / "participants_yin.csv", dtype={"eid": str})
    yang = pd.read_csv(cache / "participants_yang.csv", dtype={"eid": str})
    features = pd.read_csv(cache / "protein_features.csv")["feature"].astype(str).tolist()
    if (
        len(yin) != 37127
        or int(yin["event"].sum()) != 3442
        or len(yang) != 1766
        or len(features) != 2910
        or yin["eid"].duplicated().any()
        or yang["eid"].duplicated().any()
        or sorted(yin["outer_fold"].unique().tolist()) != [1, 2, 3, 4, 5]
    ):
        raise RuntimeError("Unified Yu participant/feature contract failed")

    def matrix(path: Path, rows: int) -> np.memmap:
        expected = rows * len(features) * 4
        if path.stat().st_size != expected:
            raise RuntimeError(f"Float32 matrix-size contract failed: {path}")
        raw = np.memmap(path, dtype="<f4", mode="r", shape=(rows * len(features),))
        return raw.reshape((rows, len(features)), order="F")

    x_yin = matrix(cache / "protein_yin.f32", len(yin))
    x_yang = matrix(cache / "protein_yang.f32", len(yang))
    yin_rows: list[pd.DataFrame] = []
    yang_rows: list[pd.DataFrame] = []
    metadata: list[dict[str, object]] = []

    for fold in range(1, 6):
        model_path = benchmark / "03_models" / f"fold{fold:02d}_yu_models" / "Yu5_ProteinAll.txt"
        if not model_path.exists():
            raise RuntimeError(f"Missing unified Yu model: {model_path}")
        booster = lgb.Booster(model_file=str(model_path))
        if len(booster.feature_name()) != len(features):
            raise RuntimeError(f"Yu fold {fold} feature count differs from locked matrix")
        test_index = np.flatnonzero(yin["outer_fold"].to_numpy() == fold)
        yin_probability = booster.predict(x_yin[test_index])
        yang_probability = booster.predict(x_yang)
        if not np.isfinite(yin_probability).all() or not np.isfinite(yang_probability).all():
            raise RuntimeError(f"Non-finite unified Yu prediction in fold {fold}")
        yin_rows.append(
            pd.DataFrame(
                {
                    "eid": yin.iloc[test_index]["eid"].to_numpy(),
                    "outer_fold": fold,
                    "time": yin.iloc[test_index]["time"].to_numpy(),
                    "event": yin.iloc[test_index]["event"].to_numpy(),
                    "series": "Yu-style LightGBM",
                    "score": logit(yin_probability),
                }
            )
        )
        yang_rows.append(
            pd.DataFrame(
                {
                    "eid": yang["eid"].to_numpy(),
                    "outer_fold": fold,
                    "series": "Yu-style LightGBM",
                    "score": logit(yang_probability),
                }
            )
        )
        metadata.append(
            {
                "model": "Yu-style LightGBM",
                "outer_fold": fold,
                "model_path": str(model_path),
                "feature_n": len(features),
                "num_trees": booster.num_trees(),
                "best_iteration": booster.best_iteration,
                "lightgbm_version": lgb.__version__,
            }
        )

    yin_output = pd.concat(yin_rows, ignore_index=True)
    yang_output = pd.concat(yang_rows, ignore_index=True)
    if (
        len(yin_output) != len(yin)
        or yin_output.duplicated(["eid"]).any()
        or len(yang_output) != 5 * len(yang)
        or yang_output.duplicated(["eid", "outer_fold"]).any()
    ):
        raise RuntimeError("Unified Yu replay output contract failed")

    Path(args.yin_output).parent.mkdir(parents=True, exist_ok=True)
    yin_output.to_csv(args.yin_output, index=False, compression="gzip")
    yang_output.to_csv(args.yang_output, index=False, compression="gzip")
    pd.DataFrame(metadata).to_csv(args.metadata_output, index=False)
    qc = {
        "status": "PASS",
        "yin_n": int(len(yin_output)),
        "yin_events": int(yin_output["event"].sum()),
        "yang_fold_rows": int(len(yang_output)),
        "feature_n": int(len(features)),
        "fold_n": 5,
    }
    Path(args.qc_output).write_text(json.dumps(qc, indent=2), encoding="utf-8")
    print(json.dumps(qc, indent=2))


if __name__ == "__main__":
    main()
