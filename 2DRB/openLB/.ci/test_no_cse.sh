#!/usr/bin/env bash

export FEATURES=DISABLE_CSE

make
make -C tests/unit build

export MAKEFLAGS="-j1"
make -C tests/unit
