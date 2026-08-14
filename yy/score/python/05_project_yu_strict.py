#!/usr/bin/env python3
"""Project a fixed native Yu LightGBM model onto the common Yin/Yang cohort."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import time
from pathlib import Path

import lightgbm as lgb
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=["preflight", "project", "status"], required=True)
    parser.add_argument("--source-root", required=True)
    parser.add_argument("--common-root", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--raw-protein-file", required=True)
    return parser.parse_args()


def canonical(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def participant_file(root: Path, side: str) -> Path:
    for name in (f"participants_{side}.csv.gz", f"participants_{side}.csv"):
        candidate = root / name
        if candidate.exists():
            return candidate
    raise RuntimeError(f"Common participant file missing for {side}: {root}")


def matrix(path: Path, rows: int, columns: int) -> np.memmap:
    expected = rows * columns * 4
    if not path.exists() or path.stat().st_size != expected:
        raise RuntimeError(f"Float32 matrix contract failed: {path}")
    raw = np.memmap(path, dtype="<f4", mode="r", shape=(rows * columns,))
    return raw.reshape((rows, columns), order="F")


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


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    args = arguments()
    source = Path(args.source_root)
    common = Path(args.common_root)
    output = Path(args.output_root)
    raw_protein_file = Path(args.raw_protein_file)
    model_file = source / "09_models" / "cad__Protein.txt"
    prediction_file = source / "09_models" / "test_predictions.csv.gz"
    yin_file = participant_file(common, "yin")
    yang_file = participant_file(common, "yang")
    feature_file = common / "protein_features.csv"
    required = [
        model_file,
        prediction_file,
        raw_protein_file,
        yin_file,
        yang_file,
        feature_file,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise RuntimeError("Yu strict projection input missing: " + "; ".join(missing))

    outputs = [
        output / "scores_yin.csv.gz",
        output / "scores_yang.csv.gz",
        output / "roc_native.csv.gz",
        output / "metrics.csv",
        output / "parameter_audit.csv",
        output / "input_manifest.csv",
        output / "COMPLETE",
    ]
    if args.stage == "status":
        print("COMPLETE" if outputs[-1].exists() else "NOT_COMPLETE")
        raise SystemExit(0 if outputs[-1].exists() else 1)

    yin = pd.read_csv(yin_file, dtype={"eid": str})
    yang = pd.read_csv(yang_file, dtype={"eid": str})
    common_features = pd.read_csv(feature_file)["feature"].astype(str).tolist()
    if (
        len(yin) != 37127
        or len(yang) != 1766
        or len(common_features) != 2910
        or yin["eid"].duplicated().any()
        or yang["eid"].duplicated().any()
    ):
        raise RuntimeError("Common Yin/Yang projection contract failed")

    booster = lgb.Booster(model_file=str(model_file))
    model_features = booster.feature_name()
    common_keys = [canonical(value) for value in common_features]
    if len(set(common_keys)) != len(common_keys):
        raise RuntimeError("Common features are ambiguous after canonicalization")
    common_lookup = {value: index for index, value in enumerate(common_keys)}
    feature_index = np.asarray([common_lookup.get(canonical(value), -1) for value in model_features], dtype=int)
    mapped_feature_index = feature_index[feature_index >= 0]
    if len(np.unique(mapped_feature_index)) != len(mapped_feature_index):
        raise RuntimeError("Yu strict feature mapping is not one-to-one")

    raw_separator = "\t" if str(raw_protein_file).lower().endswith((".tsv", ".tsv.gz", ".txt")) else ","
    raw_header = pd.read_csv(raw_protein_file, sep=raw_separator, nrows=0).columns.astype(str).tolist()
    eid_candidates = [value for value in ("eid", "id", "f.eid", "participant_id") if value in raw_header]
    if not eid_candidates:
        raise RuntimeError("Yu raw protein input has no EID column")
    raw_eid = eid_candidates[0]
    raw_keys = [canonical(value) for value in raw_header]
    duplicated_raw_keys = {value for value in raw_keys if raw_keys.count(value) > 1}
    if duplicated_raw_keys.intersection(canonical(value) for value in model_features):
        raise RuntimeError("Yu raw protein features are ambiguous after canonicalization")
    raw_lookup = {value: index for index, value in enumerate(raw_keys)}
    raw_missing = [value for value in model_features if canonical(value) not in raw_lookup]
    if raw_missing:
        raise RuntimeError("Yu strict model features missing from raw protein input: " + ", ".join(raw_missing))
    raw_features = [raw_header[raw_lookup[canonical(value)]] for value in model_features]

    print(
        f"YU_STRICT_PROJECTION_PREFLIGHT_PASS selected={len(model_features)} "
        f"common_mapped={len(mapped_feature_index)} yin={len(yin)} yang={len(yang)}"
    )
    if args.stage == "preflight":
        return
    if all(path.exists() for path in outputs):
        print("YU_STRICT_PROJECTION_PRESERVED")
        return
    if any(path.exists() for path in outputs):
        raise RuntimeError(f"Partial Yu strict projection exists; refusing overwrite: {output}")

    raw = pd.read_csv(
        raw_protein_file,
        sep=raw_separator,
        usecols=[raw_eid] + raw_features,
        dtype={raw_eid: str},
        low_memory=False,
    ).rename(columns={raw_eid: "eid"})
    raw["eid"] = raw["eid"].astype(str).str.replace(r"\.0$", "", regex=True).str.strip()
    if raw["eid"].duplicated().any():
        raise RuntimeError("Yu raw protein input contains duplicate EIDs")

    def aligned_matrix(participants: pd.DataFrame, label: str) -> pd.DataFrame:
        merged = participants[["eid"]].merge(raw, on="eid", how="left", validate="one_to_one", indicator=True)
        if merged["_merge"].ne("both").any():
            missing_eids = merged.loc[merged["_merge"].ne("both"), "eid"].astype(str).tolist()
            raise RuntimeError(f"Yu raw protein input misses {label} EIDs: {missing_eids[:10]}")
        values = merged[raw_features].copy()
        values.columns = model_features
        return values

    x_yin = aligned_matrix(yin, "Yin")
    x_yang = aligned_matrix(yang, "Yang")
    yin_probability = booster.predict(x_yin)
    yang_probability = booster.predict(x_yang)
    epsilon = 1e-7
    yin_probability = np.clip(np.asarray(yin_probability, dtype=float), epsilon, 1 - epsilon)
    yang_probability = np.clip(np.asarray(yang_probability, dtype=float), epsilon, 1 - epsilon)
    yin_raw = np.log(yin_probability / (1 - yin_probability))
    yang_raw = np.log(yang_probability / (1 - yang_probability))
    center = float(np.mean(yin_raw))
    scale = float(np.std(yin_raw, ddof=1))
    if not np.isfinite(center) or not np.isfinite(scale) or scale <= 0:
        raise RuntimeError("Yu strict projected score scale failed")

    yin_output = pd.DataFrame(
        {
            "eid": yin["eid"].astype(str),
            "method": "yu-strict",
            "cohort_side": "Yin",
            "score_raw": yin_raw,
            "score_z": (yin_raw - center) / scale,
            "score_source": "fixed native Yu protein model projected to common Yin",
        }
    )
    yang_output = pd.DataFrame(
        {
            "eid": yang["eid"].astype(str),
            "method": "yu-strict",
            "cohort_side": "Yang",
            "score_raw": yang_raw,
            "score_z": (yang_raw - center) / scale,
            "score_source": "fixed native Yu protein model projected to common Yang",
        }
    )

    native = pd.read_csv(prediction_file, dtype={"eid": str})
    native = native.loc[
        (native["outcome_id"] == "cad") & (native["model_id"] == "Protein"),
        ["eid", "y", "prediction"],
    ].copy()
    if len(native) < 1000 or native["eid"].duplicated().any() or native["y"].nunique() != 2:
        raise RuntimeError("Yu strict native prediction contract failed")
    native["eid"] = native["eid"].astype(str).str.replace(r"\.0$", "", regex=True).str.strip()
    native_matrix = aligned_matrix(native, "native held-out")
    native_replay = np.asarray(booster.predict(native_matrix), dtype=float)
    replay_difference = native_replay - native["prediction"].to_numpy(dtype=float)
    replay_correlation = float(np.corrcoef(native_replay, native["prediction"].to_numpy(dtype=float))[0, 1])
    replay_max_abs_difference = float(np.max(np.abs(replay_difference)))
    replay_rmse = float(np.sqrt(np.mean(replay_difference**2)))
    if not np.isfinite(replay_correlation) or replay_correlation < 0.999:
        raise RuntimeError(
            f"Yu frozen-score replay validation failed: correlation={replay_correlation:.8f}"
        )
    auc_value = float(roc_auc_score(native["y"], native["prediction"]))
    false_positive_rate, true_positive_rate, _ = roc_curve(native["y"], native["prediction"])
    roc_native = pd.DataFrame(
        {
            "method": "yu-strict",
            "series": "Yu strict",
            "false_positive_rate": false_positive_rate,
            "true_positive_rate": true_positive_rate,
        }
    )
    metrics = pd.DataFrame(
        {
            "method": ["yu-strict"],
            "level": ["native_heldout"],
            "n": [len(native)],
            "events": [int(native["y"].sum())],
            "auc": [auc_value],
            "auc_ci_low": [np.nan],
            "auc_ci_high": [np.nan],
            "estimand": ["native held-out eventual incident CAD binary AUC"],
            "projection_validation_n": [len(native)],
            "projection_validation_pearson": [replay_correlation],
            "projection_validation_max_abs_difference": [replay_max_abs_difference],
            "projection_validation_rmse": [replay_rmse],
        }
    )
    parameter_audit = pd.DataFrame(
        {
            "feature_order": np.arange(1, len(model_features) + 1),
            "model_feature": model_features,
            "raw_feature": raw_features,
            "common_feature": [common_features[index] if index >= 0 else "" for index in feature_index],
            "common_feature_index_1based": [index + 1 if index >= 0 else np.nan for index in feature_index],
        }
    )
    manifest = pd.DataFrame(
        {
            "role": [
                "strict_model",
                "strict_predictions",
                "raw_protein",
                "participants_yin",
                "participants_yang",
                "protein_features",
            ],
            "path": [str(path) for path in required],
            "bytes": [path.stat().st_size for path in required],
            "sha256": [file_sha256(path) for path in required],
        }
    )

    output.mkdir(parents=True, exist_ok=True)
    atomic_csv(yin_output, outputs[0])
    atomic_csv(yang_output, outputs[1])
    atomic_csv(roc_native, outputs[2])
    atomic_csv(metrics, outputs[3])
    atomic_csv(parameter_audit, outputs[4])
    atomic_csv(manifest, outputs[5])
    atomic_text(
        [
            "status=COMPLETE",
            "method=yu-strict",
            f"source_root={source}",
            f"selected_proteins={len(model_features)}",
            f"yin_n={len(yin)}",
            f"yang_n={len(yang)}",
            f"yin_center={center:.16g}",
            f"yin_scale={scale:.16g}",
            f"native_auc={auc_value:.16g}",
            f"completed_at={time.strftime('%Y-%m-%dT%H:%M:%S%z')}",
        ],
        outputs[6],
    )
    print(f"YU_STRICT_PROJECTION_COMPLETE\nOUTPUT_ROOT={output}\nNATIVE_AUC={auc_value:.6f}")


if __name__ == "__main__":
    main()
