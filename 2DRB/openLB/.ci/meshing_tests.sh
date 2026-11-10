#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
uv run script/meshing/tests/run_tests.py -q
