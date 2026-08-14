#!/usr/bin/env python3
"""Apply the frozen Yu CAD protein LightGBM model to the locked Yin cohort."""

import argparse
import json
from pathlib import Path

import lightgbm as lgb
import numpy as np
import pandas as pd


def normalize_eid(values: pd.Series) -> pd.Series:
    out = values.astype(str).str.strip()
    out = out.str.replace(r"\.0$", "", regex=True)
    return out


def load_booster(model_path: str) -> lgb.Booster:
    booster = lgb.Booster(model_file=model_path)
    if len(booster.feature_name()) != 359:
        raise RuntimeError(f"Yu feature contract failed: {len(booster.feature_name())} != 359")
    return booster


def write_features(args: argparse.Namespace) -> None:
    booster = load_booster(args.model)
    features = booster.feature_name()
    pd.DataFrame({"feature": features, "feature_order": np.arange(1, len(features) + 1)}).to_csv(
        args.features_out, index=False
    )
    print(f"Yu feature manifest written: {len(features)} proteins")


def score(args: argparse.Namespace) -> None:
    booster = load_booster(args.model)
    features = booster.feature_name()
    if bool(args.matrix) == bool(args.raw_source):
        raise RuntimeError("Provide exactly one of --matrix or --raw-source")
    if args.raw_source:
        separator = "\t" if str(args.raw_source).lower().endswith((".tsv", ".tsv.gz")) else ","
        header = pd.read_csv(args.raw_source, sep=separator, nrows=0).columns.tolist()
        lookup = {str(name).upper(): str(name) for name in header}
        if len(lookup) != len(header):
            raise RuntimeError("Raw Yu protein header is case-insensitively ambiguous")
        required = ["eid", *features]
        missing_header = [name for name in required if str(name).upper() not in lookup]
        if missing_header:
            raise RuntimeError(f"Raw Yu source lacks model columns: {missing_header[:10]}")
        source_columns = [lookup[str(name).upper()] for name in required]
        target = pd.read_csv(args.raw_source, sep=separator, usecols=source_columns)
        target = target.rename(columns={source: required[i] for i, source in enumerate(source_columns)})
    else:
        target = pd.read_csv(args.matrix)
    if "eid" not in target.columns:
        raise RuntimeError("Target Yu matrix has no eid column")
    target["eid"] = normalize_eid(target["eid"])
    if args.eid_filter:
        filter_table = pd.read_csv(args.eid_filter, compression="infer")
        if "eid" not in filter_table.columns:
            raise RuntimeError("Yu eid-filter table has no eid column")
        if {"outcome_id", "model_id"}.issubset(filter_table.columns):
            filter_table = filter_table.loc[
                (filter_table["outcome_id"].astype(str) == "cad")
                & (filter_table["model_id"].astype(str) == "Protein")
            ]
        filter_eids = set(normalize_eid(filter_table["eid"]))
        target = target.loc[target["eid"].isin(filter_eids)].copy()
        if target.shape[0] != len(filter_eids):
            raise RuntimeError(
                f"Yu raw source/filter alignment failed: rows={target.shape[0]}, requested={len(filter_eids)}"
            )
    missing = [name for name in features if name not in target.columns]
    extra = [name for name in target.columns if name != "eid" and name not in features]
    if missing or extra:
        raise RuntimeError(f"Yu matrix contract failed; missing={missing[:10]}, extra={extra[:10]}")
    if target["eid"].duplicated().any():
        raise RuntimeError("Duplicate target EIDs")

    prediction = booster.predict(target[features], num_iteration=booster.best_iteration)
    if prediction.shape[0] != target.shape[0] or not np.isfinite(prediction).all():
        raise RuntimeError("Yu target prediction contract failed")
    output = pd.DataFrame({"eid": target["eid"], "yu_probability": prediction})
    output.to_csv(args.scores_out, index=False, compression="gzip")

    if args.trusted_validation_qc:
        trusted = json.loads(Path(args.trusted_validation_qc).read_text(encoding="utf-8"))
        trusted_ok = bool(
            trusted.get("validation_pass", False)
            and int(trusted.get("feature_n", -1)) == len(features)
            and float(trusted.get("validation_pearson", 0.0)) >= 0.999999
            and float(trusted.get("validation_max_abs", float("inf"))) <= 1e-4
        )
        if not trusted_ok:
            raise RuntimeError("Trusted Yu frozen-model validation QC failed")
        qc = {
            "lightgbm_version": lgb.__version__,
            "target_n": int(output.shape[0]),
            "feature_n": int(len(features)),
            "validation_mode": "reused_locked_yin_implementation_check",
            "trusted_validation_pass": True,
            "validation_pearson": float(trusted["validation_pearson"]),
            "validation_max_abs": float(trusted["validation_max_abs"]),
            "validation_pass": True,
        }
        Path(args.qc_out).write_text(json.dumps(qc, indent=2), encoding="utf-8")
        print(json.dumps(qc, indent=2))
        return

    frozen = pd.read_csv(args.test_predictions, compression="infer")
    frozen = frozen.loc[
        (frozen["outcome_id"].astype(str) == "cad")
        & (frozen["model_id"].astype(str) == "Protein"),
        ["eid", "prediction"],
    ].copy()
    frozen["eid"] = normalize_eid(frozen["eid"])
    frozen = frozen.rename(columns={"prediction": "frozen_probability"})
    overlap = output.merge(frozen, on="eid", how="inner", validate="one_to_one")
    if overlap.shape[0] < 1000:
        raise RuntimeError(f"Yu validation overlap too small: {overlap.shape[0]}")
    delta = overlap["yu_probability"].to_numpy() - overlap["frozen_probability"].to_numpy()
    correlation = float(np.corrcoef(overlap["yu_probability"], overlap["frozen_probability"])[0, 1])
    qc = {
        "lightgbm_version": lgb.__version__,
        "target_n": int(output.shape[0]),
        "feature_n": int(len(features)),
        "validation_overlap_n": int(overlap.shape[0]),
        "validation_pearson": correlation,
        "validation_mae": float(np.mean(np.abs(delta))),
        "validation_max_abs": float(np.max(np.abs(delta))),
    }
    qc["validation_pass"] = bool(
        qc["validation_pearson"] >= 0.999999
        and qc["validation_mae"] <= 1e-6
        and qc["validation_max_abs"] <= 1e-4
    )
    Path(args.qc_out).write_text(json.dumps(qc, indent=2), encoding="utf-8")
    print(json.dumps(qc, indent=2))
    if not qc["validation_pass"]:
        raise RuntimeError("Yu frozen-model validation failed")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=["features", "score"], required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--features-out")
    parser.add_argument("--matrix")
    parser.add_argument("--raw-source")
    parser.add_argument("--eid-filter")
    parser.add_argument("--scores-out")
    parser.add_argument("--test-predictions")
    parser.add_argument("--trusted-validation-qc")
    parser.add_argument("--qc-out")
    args = parser.parse_args()
    if args.stage == "features":
        if not args.features_out:
            parser.error("--features-out is required for features stage")
        write_features(args)
    else:
        required = [args.scores_out, args.qc_out]
        if any(value is None for value in required):
            parser.error("--scores-out and --qc-out are required")
        if bool(args.matrix) == bool(args.raw_source):
            parser.error("Provide exactly one of --matrix or --raw-source")
        if bool(args.test_predictions) == bool(args.trusted_validation_qc):
            parser.error("Provide exactly one of --test-predictions or --trusted-validation-qc")
        score(args)


if __name__ == "__main__":
    main()
