#!/usr/bin/env python
"""Source-locked Yu/Chen 2025 LightGBM selection and prediction stages."""

import argparse
import hashlib
import json
import os
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMClassifier
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import roc_auc_score, roc_curve


def sha_text(values):
    return hashlib.sha256("\n".join(map(str, values)).encode("utf-8")).hexdigest()


def read_json(path):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def read_table(path, usecols=None):
    sep = "\t" if str(path).lower().endswith((".tsv", ".txt", ".tsv.gz")) else ","
    return pd.read_csv(path, sep=sep, usecols=usecols, low_memory=False)


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def unique_in_order(values):
    seen = set()
    ordered = []
    for value in values:
        item = str(value).strip()
        if item and item not in seen:
            seen.add(item)
            ordered.append(item)
    return ordered


def resolve_custom_panel(raw_protein_file, panel_file="", inline_proteins="", mapping_file=""):
    values = []
    source = {"file": "", "file_sha256": "", "file_column": "", "inline_n": 0}
    if panel_file:
        path = Path(panel_file)
        if not path.exists():
            raise RuntimeError(f"Custom protein panel file not found: {path}")
        table = read_table(path)
        aliases = ("feature_id", "local_feature", "protein", "assay")
        column = next((name for name in aliases if name in table.columns), None)
        if column is None:
            if table.shape[1] != 1:
                raise RuntimeError(
                    "Custom protein panel file must contain feature_id, local_feature, "
                    "protein, assay, or exactly one column"
                )
            column = table.columns[0]
        values.extend(table[column].dropna().astype(str).tolist())
        source.update({
            "file": str(path.resolve()),
            "file_sha256": file_sha256(path),
            "file_column": str(column),
        })
    inline = [item.strip() for item in str(inline_proteins).split(",") if item.strip()]
    values.extend(inline)
    source["inline_n"] = len(inline)
    requested = unique_in_order(values)
    if not requested:
        raise RuntimeError("Custom protein panel is empty")

    raw_sep = "\t" if str(raw_protein_file).lower().endswith((".tsv", ".txt", ".tsv.gz")) else ","
    raw_header = list(pd.read_csv(raw_protein_file, sep=raw_sep, nrows=0).columns)
    raw_by_lower = {}
    for feature in raw_header:
        raw_by_lower.setdefault(str(feature).lower(), []).append(str(feature))

    mapping = None
    if mapping_file:
        mapping_path = Path(mapping_file)
        if not mapping_path.exists():
            raise RuntimeError(f"Protein mapping file not found: {mapping_path}")
        mapping = read_table(mapping_path).fillna("")

    resolved = []
    request_map = []
    missing = []
    ambiguous = []
    mapping_match_columns = ("HGNC.symbol", "Assay", "OlinkID", "UniProt", "UKBPPP_ProteinID")
    for item in requested:
        direct = raw_by_lower.get(item.lower(), [])
        candidates = unique_in_order(direct)
        method = "raw_feature_case_insensitive"
        if not candidates and mapping is not None:
            matched = pd.Series(False, index=mapping.index)
            for column in mapping_match_columns:
                if column in mapping.columns:
                    matched |= mapping[column].astype(str).str.lower().eq(item.lower())
            candidate_rows = mapping.loc[matched]
            mapped = []
            if "Assay" in candidate_rows.columns:
                for assay in candidate_rows["Assay"].astype(str):
                    mapped.extend(raw_by_lower.get(assay.lower(), []))
            candidates = unique_in_order(mapped)
            method = "panel_mapping_unique_assay"
        if len(candidates) == 1:
            resolved.append(candidates[0])
            request_map.append({"requested": item, "feature_id": candidates[0], "method": method})
        elif len(candidates) > 1:
            ambiguous.append({"requested": item, "candidates": candidates})
        else:
            missing.append(item)
    if ambiguous:
        raise RuntimeError(
            "Custom protein identifiers map to multiple assays; use an exact local feature ID. "
            f"First={ambiguous[:5]}"
        )
    if missing:
        raise RuntimeError(
            "Custom protein panel contains identifiers absent from the raw protein table/mapping; "
            f"missing={len(missing)}, first={missing[:10]}. Use a local feature ID or uniquely mapped identifier."
        )
    features = unique_in_order(resolved)
    source["feature_n"] = len(features)
    source["feature_hash"] = sha_text(features)
    source["requested_n"] = len(requested)
    source["requested_to_feature"] = request_map
    source["mapping_file"] = str(Path(mapping_file).resolve()) if mapping_file else ""
    return features, source


def normalize_eid(values):
    return values.astype(str).str.replace(r"\.0$", "", regex=True).str.strip()


def build_model(params, workers, seed):
    return LGBMClassifier(
        objective="binary",
        n_estimators=int(params["n_estimators"]),
        max_depth=int(params["max_depth"]),
        num_leaves=int(params["num_leaves"]),
        subsample=float(params["subsample"]),
        subsample_freq=int(params.get("subsample_freq", 1)),
        learning_rate=float(params["learning_rate"]),
        colsample_bytree=float(params["colsample_bytree"]),
        random_state=int(seed),
        n_jobs=int(workers),
        verbosity=-1,
    )


def youden_threshold(y, prediction):
    fpr, tpr, thresholds = roc_curve(y, prediction)
    valid = np.isfinite(thresholds)
    if not valid.any():
        return 0.5
    idx = np.argmax((tpr - fpr)[valid])
    return float(thresholds[valid][idx])


def endpoint_ids(project_dir, endpoint_subset):
    outcomes = pd.read_csv(Path(project_dir) / "f" / "config" / "outcomes.csv")
    if endpoint_subset.lower() == "all":
        return outcomes.outcome_id.tolist()
    requested = [x.strip() for x in endpoint_subset.split(",") if x.strip()]
    missing = sorted(set(requested) - set(outcomes.outcome_id))
    if missing:
        raise RuntimeError(f"Unknown endpoint IDs: {missing}")
    return requested


def load_cohorts(analysis_dir):
    cohort_dir = Path(analysis_dir) / "04_cohort"
    derivation = pd.read_csv(cohort_dir / "derivation_cohort.csv.gz", dtype={"eid": str})
    test = pd.read_csv(cohort_dir / "test_cohort.csv.gz", dtype={"eid": str})
    folds = pd.read_csv(cohort_dir / "foldid.csv", dtype={"eid": str})
    derivation["eid"] = normalize_eid(derivation.eid)
    test["eid"] = normalize_eid(test.eid)
    folds["eid"] = normalize_eid(folds.eid)
    derivation = derivation.merge(folds, on="eid", how="left", validate="one_to_one")
    if derivation.foldid.isna().any():
        raise RuntimeError("Missing derivation fold IDs")
    return derivation, test


def load_proteins(raw_file, feature_ids, cohort):
    sample = pd.read_csv(raw_file, sep="\t" if str(raw_file).lower().endswith((".tsv", ".tsv.gz", ".txt")) else ",", nrows=0)
    eid_candidates = [x for x in ("eid", "id", "f.eid", "participant_id") if x in sample.columns]
    if not eid_candidates:
        raise RuntimeError("No EID column in raw protein table")
    eid_col = eid_candidates[0]
    missing = sorted(set(feature_ids) - set(sample.columns))
    if missing:
        raise RuntimeError(f"Raw protein table is missing {len(missing)} requested features; first={missing[:10]}")
    protein = read_table(raw_file, usecols=[eid_col] + list(feature_ids))
    protein = protein.rename(columns={eid_col: "eid"})
    protein["eid"] = normalize_eid(protein.eid)
    if protein.eid.duplicated().any():
        raise RuntimeError("Duplicate EIDs in raw protein table")
    cohort_x = cohort[["eid"]].merge(
        protein, on="eid", how="left", validate="one_to_one", indicator="__protein_merge"
    )
    absent = cohort_x["__protein_merge"].ne("both")
    if absent.any():
        missing_eids = cohort_x.loc[absent, "eid"].astype(str).tolist()
        raise RuntimeError(
            "Raw protein input has no row for "
            f"{len(missing_eids)} requested cohort EIDs; first={missing_eids[:10]}"
        )
    cohort_x = cohort_x.drop(columns="__protein_merge")
    missing = cohort_x[feature_ids].isna()
    cohort_x.attrs["protein_coverage"] = {
        "participant_n": int(len(cohort_x)),
        "requested_feature_n": int(len(feature_ids)),
        "all_feature_missing_n": int(missing.all(axis=1).sum()),
        "any_feature_missing_n": int(missing.any(axis=1).sum()),
        "cell_missing_fraction": float(missing.to_numpy().mean()),
        "missing_value_policy": "LightGBM native missing-value handling",
    }
    return cohort_x


def select_stage(args, cfg):
    analysis = Path(args.analysis_dir)
    selection_dir = analysis / "08_selection"
    selection_dir.mkdir(parents=True, exist_ok=True)
    candidate_file = selection_dir / "derivation_bonferroni_candidate_union.csv"
    association_file = analysis / "05_cox" / "derivation_associations.csv.gz"
    if not candidate_file.exists() or not association_file.exists():
        raise RuntimeError("Run the R cox stage first")
    candidates = pd.read_csv(candidate_file).feature_id.dropna().drop_duplicates().tolist()
    if not candidates:
        raise RuntimeError("No derivation Bonferroni candidates")
    derivation, _ = load_cohorts(analysis)
    derivation_x = load_proteins(args.raw_protein_file, candidates, derivation)
    candidate_coverage = derivation_x.attrs["protein_coverage"]
    all_candidate_missing = derivation_x[candidates].isna().all(axis=1)
    candidate_missing_eids = derivation.loc[
        all_candidate_missing,
        ["eid"] + [column for column in derivation.columns if column.startswith("event_")],
    ].copy()
    candidate_missing_eids.to_csv(
        selection_dir / "candidate_panel_all_missing_participants.csv", index=False
    )
    with open(selection_dir / "candidate_panel_missingness_qc.json", "w", encoding="utf-8") as handle:
        json.dump({"status": "PASS", **candidate_coverage}, handle, indent=2)
    endpoints = endpoint_ids(args.project_dir, args.endpoint_subset)
    params = cfg["lightgbm"]
    fraction = float(cfg["importance_cumulative_fraction"])

    importance_rows = []
    selected_rows = []
    for offset, endpoint in enumerate(endpoints):
        y = derivation[f"event_{endpoint}"].astype(int).to_numpy()
        if np.unique(y).size < 2:
            raise RuntimeError(f"Endpoint {endpoint} has only one class in derivation")
        model = build_model(params, args.workers, int(params["random_state"]) + offset)
        model.fit(derivation_x[candidates], y)
        gain = model.booster_.feature_importance(importance_type="gain").astype(float)
        total = gain.sum()
        if not np.isfinite(total) or total <= 0:
            raise RuntimeError(f"Endpoint {endpoint} produced zero total information gain")
        table = pd.DataFrame({"outcome_id": endpoint, "feature_id": candidates, "gain": gain})
        table["normalized_gain"] = table.gain / total
        table = table.sort_values(["normalized_gain", "feature_id"], ascending=[False, True]).reset_index(drop=True)
        table["rank"] = np.arange(1, len(table) + 1)
        table["cumulative_gain"] = table.normalized_gain.cumsum()
        crossing = int(np.searchsorted(table.cumulative_gain.to_numpy(), fraction, side="left"))
        table["selected_to_30pct"] = table.index <= crossing
        importance_rows.append(table)
        selected_rows.append(table.loc[table.selected_to_30pct, ["outcome_id", "feature_id", "rank", "normalized_gain", "cumulative_gain"]])

    importance = pd.concat(importance_rows, ignore_index=True)
    selected = pd.concat(selected_rows, ignore_index=True)
    local_union = sorted(selected.feature_id.unique())
    official_file = Path(args.project_dir) / "f" / "config" / "yu_cad_257_official.csv"
    if not official_file.exists():
        raise RuntimeError(f"Published S12 panel mapping is missing: {official_file}")
    official = pd.read_csv(official_file)
    required_official_columns = {"protein", "rank", "local_feature"}
    if not required_official_columns.issubset(official.columns):
        raise RuntimeError(f"Published panel mapping lacks columns: {sorted(required_official_columns - set(official.columns))}")
    official = official.sort_values("rank").drop_duplicates("local_feature")
    published_union = official.local_feature.dropna().astype(str).tolist()
    if len(published_union) != 257 or official.protein.nunique() != 257:
        raise RuntimeError(
            f"Published S12 panel contract failed: features={len(published_union)} proteins={official.protein.nunique()}"
        )
    raw_sep = "\t" if str(args.raw_protein_file).lower().endswith((".tsv", ".txt", ".tsv.gz")) else ","
    raw_header = set(pd.read_csv(args.raw_protein_file, sep=raw_sep, nrows=0).columns)
    missing_published = sorted(set(published_union) - raw_header)
    if missing_published:
        raise RuntimeError(
            f"Raw protein table is missing {len(missing_published)} published S12 features; first={missing_published[:10]}"
        )
    custom_union = []
    custom_source = {}
    if args.prediction_panel_mode == "published_257":
        union = published_union
        final_source = "official_supplementary_table_S12"
    elif args.prediction_panel_mode == "custom":
        custom_union, custom_source = resolve_custom_panel(
            args.raw_protein_file, args.custom_protein_panel_file, args.custom_proteins,
            args.panel_mapping_file
        )
        union = custom_union
        final_source = "user_custom_panel"
    else:
        union = local_union
        final_source = "local_derivation_reselection"
    importance.to_csv(selection_dir / "preliminary_lgbm_importance.csv.gz", index=False)
    selected.to_csv(selection_dir / "endpoint_selected_to_30pct.csv", index=False)
    pd.DataFrame({"feature_id": local_union}).to_csv(
        selection_dir / "local_reselected_cross_endpoint_protein_union.csv", index=False
    )
    official.to_csv(selection_dir / "published_257_panel_mapping.csv", index=False)
    pd.DataFrame({"feature_id": published_union}).to_csv(
        selection_dir / "published_257_cross_endpoint_protein_union.csv", index=False
    )
    if custom_union:
        pd.DataFrame({"feature_id": custom_union}).to_csv(
            selection_dir / "custom_cross_endpoint_protein_union.csv", index=False
        )
        pd.DataFrame(custom_source["requested_to_feature"]).to_csv(
            selection_dir / "custom_protein_identifier_mapping.csv", index=False
        )
    pd.DataFrame({"feature_id": union}).to_csv(selection_dir / "final_cross_endpoint_protein_union.csv", index=False)
    overlap_n = len(set(local_union) & set(published_union))
    with open(selection_dir / "selection_summary.json", "w", encoding="utf-8") as handle:
        json.dump({
            "status": "PASS",
            "endpoint_ids": endpoints,
            "endpoint_count": len(endpoints),
            "candidate_union_n": len(candidates),
            "final_union_n": len(union),
            "prediction_panel_mode": args.prediction_panel_mode,
            "final_panel_source": final_source,
            "local_reselected_union_n": len(local_union),
            "published_union_n": len(published_union),
            "local_published_overlap_n": overlap_n,
            "local_published_jaccard": overlap_n / len(set(local_union) | set(published_union)),
            "candidate_hash": sha_text(candidates),
            "final_union_hash": sha_text(union),
            "local_reselected_union_hash": sha_text(local_union),
            "published_union_hash": sha_text(published_union),
            "custom_panel": custom_source,
            "cumulative_gain_fraction": fraction,
            "published_anchors": {"candidate_union_n": 671, "final_union_n": 257},
            "interpretation": {
                "local_reselected": (
                    "The prediction benchmark uses derivation-only Yu-style local reselection."
                ),
                "published_257": (
                    "The prediction benchmark uses the exact published S12 panel as a reference-mode analysis."
                ),
                "custom": (
                    "The prediction benchmark uses a user-specified panel of exact local protein feature IDs."
                ),
            }[args.prediction_panel_mode],
        }, handle, indent=2)
    print(json.dumps({
        "status": "PASS", "mode": "select", "candidates": len(candidates),
        "local_reselected_union": len(local_union), "final_union": len(union),
        "prediction_panel_mode": args.prediction_panel_mode,
    }))


def train_stage(args, cfg):
    analysis = Path(args.analysis_dir)
    selection_dir = analysis / "08_selection"
    model_dir = analysis / "09_models"
    model_dir.mkdir(parents=True, exist_ok=True)
    union_file = selection_dir / "final_cross_endpoint_protein_union.csv"
    if not union_file.exists():
        raise RuntimeError("Run select first")
    features = pd.read_csv(union_file).feature_id.dropna().drop_duplicates().tolist()
    derivation, test = load_cohorts(analysis)
    derivation_x = load_proteins(args.raw_protein_file, features, derivation)
    test_x = load_proteins(args.raw_protein_file, features, test)
    panel_coverage = []
    all_missing_tables = []
    for split_name, cohort_meta, protein_frame in (
        ("derivation", derivation, derivation_x),
        ("test", test, test_x),
    ):
        coverage = dict(protein_frame.attrs["protein_coverage"])
        coverage["split"] = split_name
        panel_coverage.append(coverage)
        all_missing = protein_frame[features].isna().all(axis=1)
        if all_missing.any():
            missing_rows = cohort_meta.loc[
                all_missing,
                ["eid"] + [column for column in cohort_meta.columns if column.startswith("event_")],
            ].copy()
            missing_rows.insert(1, "split", split_name)
            all_missing_tables.append(missing_rows)
    pd.DataFrame(panel_coverage).to_csv(model_dir / "model_panel_missingness_qc.csv", index=False)
    if all_missing_tables:
        pd.concat(all_missing_tables, ignore_index=True).to_csv(
            model_dir / "model_panel_all_missing_participants.csv", index=False
        )
    else:
        pd.DataFrame(columns=["eid", "split"]).to_csv(
            model_dir / "model_panel_all_missing_participants.csv", index=False
        )
    derivation_x = derivation[["eid", "foldid", "score2_raw"] + [c for c in derivation.columns if c.startswith("event_")]].merge(
        derivation_x, on="eid", how="left", validate="one_to_one"
    )
    test_x = test[["eid", "score2_raw"] + [c for c in test.columns if c.startswith("event_")]].merge(
        test_x, on="eid", how="left", validate="one_to_one"
    )
    endpoints = endpoint_ids(args.project_dir, args.endpoint_subset)
    for stale_model in model_dir.glob("*__*.txt"):
        stale_model.unlink()
    for stage in ("evaluate", "figures", "report"):
        marker = analysis / "00_logs" / f"{stage}.done.json"
        if marker.exists():
            marker.unlink()
    params = cfg["lightgbm"]
    prediction_rows = []
    importance_rows = []
    fold_rows = []
    model_rows = []
    design_rows = []
    model_jobs = max(1, min(int(args.model_jobs), 5))
    workers_per_model = max(1, int(args.workers) // model_jobs)

    specs = {
        "SCORE2": ["SCORE2_calibrated"],
        "Protein": features,
        "Protein_SCORE2": features + ["SCORE2_calibrated"],
    }

    seed_group = {
        "SCORE2": 0,
        "Protein": 1,
        "Protein_SCORE2": 2,
    }

    for endpoint_offset, endpoint in enumerate(endpoints):
        y = derivation_x[f"event_{endpoint}"].astype(int).to_numpy()
        y_test = test_x[f"event_{endpoint}"].astype(int).to_numpy()
        oof_score = np.full(len(derivation_x), np.nan)
        oof_models = {name: np.full(len(derivation_x), np.nan) for name in specs}
        for fold in sorted(derivation_x.foldid.unique()):
            tr = derivation_x.foldid.to_numpy() != fold
            va = ~tr
            iso = IsotonicRegression(out_of_bounds="clip").fit(derivation_x.loc[tr, "score2_raw"], y[tr])
            score_train = iso.transform(derivation_x.loc[tr, "score2_raw"])
            score_valid = iso.transform(derivation_x.loc[va, "score2_raw"])
            train_fold = derivation_x.loc[tr].copy()
            valid_fold = derivation_x.loc[va].copy()
            train_fold["SCORE2_calibrated"] = score_train
            valid_fold["SCORE2_calibrated"] = score_valid
            oof_score[va] = score_valid
            def fit_fold_model(item):
                model_name, columns = item
                seed = int(params["random_state"]) + endpoint_offset * 100 + int(fold) * 10 + seed_group[model_name]
                model = build_model(params, workers_per_model, seed)
                model.fit(train_fold[columns], y[tr])
                pred = model.predict_proba(valid_fold[columns])[:, 1]
                return model_name, pred

            with ThreadPoolExecutor(max_workers=model_jobs) as executor:
                fold_results = list(executor.map(fit_fold_model, specs.items()))
            for model_name, pred in fold_results:
                oof_models[model_name][va] = pred
                fold_rows.append({
                    "outcome_id": endpoint, "model_id": model_name, "fold": int(fold),
                    "n": int(va.sum()), "events": int(y[va].sum()), "auc": float(roc_auc_score(y[va], pred)),
                })

        final_iso = IsotonicRegression(out_of_bounds="clip").fit(derivation_x.score2_raw, y)
        derivation_x["SCORE2_calibrated"] = final_iso.transform(derivation_x.score2_raw)
        test_x["SCORE2_calibrated"] = final_iso.transform(test_x.score2_raw)
        def fit_final_model(item):
            model_name, columns = item
            threshold = youden_threshold(y, oof_models[model_name])
            seed = int(params["random_state"]) + endpoint_offset * 100 + seed_group[model_name]
            model = build_model(params, workers_per_model, seed)
            model.fit(derivation_x[columns], y)
            prediction = model.predict_proba(test_x[columns])[:, 1]
            gain = model.booster_.feature_importance(importance_type="gain").astype(float)
            return model_name, columns, threshold, model, prediction, gain

        with ThreadPoolExecutor(max_workers=model_jobs) as executor:
            final_results = list(executor.map(fit_final_model, specs.items()))
        for model_name, columns, threshold, model, prediction, gain in final_results:
            safe_name = f"{endpoint}__{model_name}"
            model.booster_.save_model(model_dir / f"{safe_name}.txt")
            prediction_rows.append(pd.DataFrame({
                "eid": test_x.eid, "outcome_id": endpoint, "y": y_test,
                "model_id": model_name, "prediction": prediction, "threshold": threshold,
            }))
            importance_rows.append(pd.DataFrame({
                "outcome_id": endpoint, "model_id": model_name, "feature": columns,
                "gain": gain, "standardized_gain": gain / gain.sum() if gain.sum() > 0 else np.nan,
            }))
            model_rows.append({
                "outcome_id": endpoint, "model_id": model_name, "feature_n": len(columns),
                "oof_auc": float(roc_auc_score(y, oof_models[model_name])),
                "threshold_source": "derivation_oof_youden", "threshold": threshold,
            })
            design_rows.append({
                "outcome_id": endpoint,
                "model_id": model_name,
                "feature_n": len(columns),
                "feature_hash": sha_text(columns),
            })

    pd.concat(prediction_rows, ignore_index=True).to_csv(model_dir / "test_predictions.csv.gz", index=False)
    pd.concat(importance_rows, ignore_index=True).to_csv(model_dir / "final_model_importance.csv.gz", index=False)
    pd.DataFrame(fold_rows).to_csv(model_dir / "derivation_fold_metrics.csv", index=False)
    pd.DataFrame(model_rows).to_csv(model_dir / "model_summary.csv", index=False)
    pd.DataFrame(design_rows).to_csv(model_dir / "model_design_contract.csv", index=False)
    with open(model_dir / "training_manifest.json", "w", encoding="utf-8") as handle:
        json.dump({
            "status": "PASS", "endpoint_ids": endpoints,
            "endpoint_count": len(endpoints), "protein_n": len(features),
            "protein_hash": sha_text(features), "parameters": params,
            "prediction_panel_mode": args.prediction_panel_mode,
            "score2_isotonic_fit_in_derivation_only": True,
            "test_predictions_not_used_for_selection_or_tuning": True,
            "model_jobs": model_jobs,
            "workers_per_model": workers_per_model,
            "lightgbm_version": __import__("lightgbm").__version__,
            "sklearn_version": __import__("sklearn").__version__,
        }, handle, indent=2)
    print(json.dumps({
        "status": "PASS", "mode": "train", "endpoints": len(endpoints),
        "protein_n": len(features),
    }))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["select", "train"], required=True)
    parser.add_argument("--project-dir", required=True)
    parser.add_argument("--analysis-dir", required=True)
    parser.add_argument("--raw-protein-file", required=True)
    parser.add_argument("--panel-mapping-file", default="")
    parser.add_argument("--config", required=True)
    parser.add_argument("--endpoint-subset", default="all")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--model-jobs", type=int, default=3)
    parser.add_argument(
        "--prediction-panel-mode", choices=["published_257", "local_reselected", "custom"],
        default="published_257"
    )
    parser.add_argument("--custom-protein-panel-file", default="")
    parser.add_argument("--custom-proteins", default="")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    cfg = read_json(args.config)
    required_python = str(cfg.get("required_python_major_minor", "3.9"))
    actual_python = f"{sys.version_info.major}.{sys.version_info.minor}"
    if actual_python != required_python:
        raise RuntimeError(f"Formal reproduction requires Python {required_python}; found {actual_python}")
    required_lightgbm = str(cfg.get("required_lightgbm_version", "3.3.2"))
    actual_lightgbm = __import__("lightgbm").__version__
    if actual_lightgbm != required_lightgbm:
        raise RuntimeError(f"Formal reproduction requires lightgbm {required_lightgbm}; found {actual_lightgbm}")
    summary_path = Path(args.analysis_dir) / ("08_selection/selection_summary.json" if args.mode == "select" else "09_models/training_manifest.json")
    endpoints = endpoint_ids(args.project_dir, args.endpoint_subset)
    if args.prediction_panel_mode == "custom":
        expected_custom, _ = resolve_custom_panel(
            args.raw_protein_file, args.custom_protein_panel_file, args.custom_proteins,
            args.panel_mapping_file
        )
        expected_custom_hash = sha_text(expected_custom)
    else:
        expected_custom_hash = None
    if args.resume and summary_path.exists():
        summary = read_json(summary_path)
        endpoint_matches = summary.get("endpoint_ids") == endpoints
        mode_matches = summary.get("prediction_panel_mode") == args.prediction_panel_mode and endpoint_matches
        if args.mode == "select":
            union_file = Path(args.analysis_dir) / "08_selection/final_cross_endpoint_protein_union.csv"
            output_matches = union_file.exists()
            if args.prediction_panel_mode == "custom":
                output_matches = output_matches and summary.get("final_union_hash") == expected_custom_hash
        else:
            union_file = Path(args.analysis_dir) / "08_selection/final_cross_endpoint_protein_union.csv"
            prediction_file = Path(args.analysis_dir) / "09_models/test_predictions.csv.gz"
            output_matches = union_file.exists() and prediction_file.exists()
            if output_matches:
                current_features = pd.read_csv(union_file).feature_id.dropna().drop_duplicates().tolist()
                output_matches = summary.get("protein_hash") == sha_text(current_features)
        if summary.get("status") == "PASS" and mode_matches and output_matches:
            print(json.dumps({"status": "PASS", "mode": args.mode, "resume": "skipped_validated_stage"}))
            return
    if args.mode == "select":
        select_stage(args, cfg)
    else:
        train_stage(args, cfg)


if __name__ == "__main__":
    main()
