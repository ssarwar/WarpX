#!/usr/bin/env python3
"""Record source hashes, including new fixtures and literal symlink targets."""

import hashlib
import json
import subprocess
import sys
from pathlib import Path

paths = (
    subprocess.check_output(
        [
            "git",
            "ls-files",
            "-z",
            "--cached",
            "--others",
            "--exclude-standard",
            "--",
            "Source",
            "Python",
            "Examples/Tests",
            "Examples/analysis_default_restart.py",
            "tests/unit",
            "Tools/Algorithms",
            "Tools/machines",
            "Docs",
            "CMakeLists.txt",
            "cmake",
            "dependencies.json",
            "requirements.txt",
        ]
    )
    .decode()
    .split("\0")
)
manifest = {}
for name in sorted(set(paths) - {""}):
    path = Path(name)
    if path.is_symlink():
        target = str(path.readlink())
        manifest[name] = dict(
            kind="symlink",
            target=target,
            sha256=hashlib.sha256(target.encode()).hexdigest(),
        )
    else:
        manifest[name] = dict(
            kind="file", sha256=hashlib.sha256(path.read_bytes()).hexdigest()
        )
Path(sys.argv[1]).write_text(json.dumps(manifest, indent=2) + "\n")
