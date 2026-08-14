import ast
import importlib.util
import subprocess
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[2]
RUNNER = PROJECT / "f" / "tools" / "yu_runner.py"

shell_files = list(PROJECT.rglob("*.sh"))
assert shell_files == [PROJECT / "yu.sh"], shell_files
assert not list(PROJECT.rglob("*.ps1"))

ast.parse(RUNNER.read_text(encoding="utf-8"))
spec = importlib.util.spec_from_file_location("yu_runner", RUNNER)
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)
assert runner.valid_model_runtime("3.9.23 3.3.2\n")
assert not runner.valid_model_runtime("3.10.1 3.3.2\n")
subprocess.run(["bash", "-n", str(PROJECT / "yu.sh")], check=True)
help_text = subprocess.check_output(
    ["bash", str(PROJECT / "yu.sh"), "--help"], text=True
)
for token in (
    "setup", "doctor", "status", "package", "finalize", "1-4",
    "--analysis-project", "--genotype-root", "--model-python",
):
    assert token in help_text, token

bad = subprocess.run(
    ["bash", str(PROJECT / "yu.sh"), "not-a-command"],
    text=True, capture_output=True,
)
assert bad.returncode != 0

print("SINGLE ENTRYPOINT TESTS PASSED")
