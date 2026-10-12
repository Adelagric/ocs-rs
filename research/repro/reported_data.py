"""Validate archived summary hashes before explicitly requested reported plots."""
import hashlib
import json
from pathlib import Path


def reported_directory():
    directory = Path(__file__).resolve().parent / "reported/v0.4.0"
    manifest = json.loads((directory / "manifest.json").read_text())
    for name, expected in manifest["sha256"].items():
        if hashlib.sha256((directory / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"historical summary changed: {name}")
    return directory, manifest
