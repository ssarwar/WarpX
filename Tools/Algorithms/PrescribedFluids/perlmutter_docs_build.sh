#!/usr/bin/env bash
# Render the complete documentation in an isolated Perlmutter environment.
set -eo pipefail
root=$(git rev-parse --show-toplevel)
audit="$root/build/reaudit-2026-09-21"
source "$audit/perlmutter_gpu_warpx.profile"
venv="$audit/software/venvs/warpx-docs"
if [[ ! -f "$venv/bin/activate" ]]; then
    python3 -m venv "$venv"
fi
source "$venv/bin/activate"
# Import the pure Python interface directly; keep the CUDA runtime environment
# independent from Sphinx's pinned dependencies.
sed '/^-e /d' "$root/Docs/requirements.txt" > "$audit/docs-requirements.txt"
python3 -m pip install --timeout 60 -r "$root/requirements.txt" \
    -r "$audit/docs-requirements.txt" > "$audit/docs-install.log" 2>&1
export PYTHONPATH="$root/Python"
python3 -m sphinx -E -b html "$root/Docs/source" "$audit/docs-html" \
    > "$audit/docs-build.log" 2>&1
python3 -m pip freeze > "$audit/docs-packages.txt"
