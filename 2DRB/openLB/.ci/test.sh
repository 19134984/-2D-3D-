#!/usr/bin/env bash

set -euo pipefail

make
make -C tests/unit build
export MAKEFLAGS="-j1"
make -C tests/unit
