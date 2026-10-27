"""Read-only exact replay; identities and algorithm tests remain separate from proof."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def main():
    root = Path(__file__).resolve().parent
    manifest = json.loads((root / "manifest.json").read_text())
    expected = {"README.md", "sources.md", "proof.md", "status.json", "LICENSE", "requirements.txt", "witness.json", "check.py", "gamma_cumulants.py", "test_gamma_cumulants.py", "reproduce.py"}
    if set(manifest["files"]) != expected:
        raise ValueError("manifest file-set mismatch")
    for name, digest in manifest["files"].items():
        if hashlib.sha256((root / name).read_bytes()).hexdigest() != digest:
            raise ValueError("integrity mismatch: " + name)
    for name in ["check.py", "test_gamma_cumulants.py"]:
        subprocess.run([sys.executable, "-E", "-B", str(root / name)], check=True)
    print("Exact cumulant packet replay passed; full-domain theorem rests on proof.md.")


if __name__ == "__main__":
    main()
