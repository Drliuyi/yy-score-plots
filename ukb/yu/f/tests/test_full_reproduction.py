import ast
import hashlib
import tempfile
from pathlib import Path

import pandas as pd


F_DIR = Path(__file__).resolve().parents[1]
PROJECT = F_DIR.parent
MODEL = F_DIR / "python" / "04_full_reproduction.py"
RUNNER = F_DIR / "tools" / "yu_runner.py"

model_text = MODEL.read_text(encoding="utf-8")
runner_text = RUNNER.read_text(encoding="utf-8")
model_tree = ast.parse(model_text)
ast.parse(runner_text)

for token in (
    "derivation_bonferroni_candidate_union.csv",
    "selected_to_30pct",
    "final_cross_endpoint_protein_union.csv",
    "SCORE2",
    "Protein",
    "Protein_SCORE2",
    "ThreadPoolExecutor",
    "test_predictions_not_used_for_selection_or_tuning",
    "LightGBM native missing-value handling",
    'choices=["published_257", "local_reselected", "custom"]',
):
    assert token in model_text, token

for retired in ("YYScore", "yys_mode", "--yys-mode"):
    assert retired not in model_text, retired

for token in (
    '"core": list(range(1, 5))',
    '"downstream": list(range(5, 11))',
    '"all": list(range(1, 12))',
    '"figures": [11]',
    '"finalize": [10, 11]',
    '"Z:/projects/genotype_pc_nas/imputed_pgen_autosomes"',
    '"yu_proteomic_repo_v3 is protected',
    'PRS_SCORER',
    'LightGBM 3.3.2',
):
    assert token in runner_text, token

nodes = {
    node.name: node for node in model_tree.body
    if isinstance(node, ast.FunctionDef) and node.name in {
        "sha_text", "read_table", "file_sha256", "unique_in_order",
        "resolve_custom_panel", "normalize_eid", "load_proteins",
    }
}
namespace = {"Path": Path, "pd": pd, "hashlib": hashlib}
exec(compile(ast.Module(body=list(nodes.values()), type_ignores=[]), str(MODEL), "exec"), namespace)

with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    raw = root / "protein.tsv"
    panel = root / "panel.csv"
    pd.DataFrame({"eid": ["1", "2"], "gdf15": [1.0, 2.0], "nppb": [2.0, 3.0]}).to_csv(
        raw, sep="\t", index=False
    )
    pd.DataFrame({"feature_id": ["GDF15", "NPPB", "GDF15"]}).to_csv(panel, index=False)
    features, manifest = namespace["resolve_custom_panel"](raw, panel)
    assert features == ["gdf15", "nppb"]
    assert manifest["feature_n"] == 2

    cohort = pd.DataFrame({"eid": ["1", "2"]})
    loaded = namespace["load_proteins"](raw, features, cohort)
    assert len(loaded) == 2
    assert loaded.attrs["protein_coverage"]["all_feature_missing_n"] == 0

    try:
        namespace["load_proteins"](raw, features, pd.DataFrame({"eid": ["3"]}))
    except RuntimeError as error:
        assert "has no row" in str(error)
    else:
        raise AssertionError("Absent protein EID did not fail")

print("FULL REPRODUCTION STATIC TESTS PASSED")
