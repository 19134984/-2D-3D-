#!/usr/bin/env bash

make

pushd tests/integration/cavity2d
make
./cavity2d 100
popd

pushd tests/integration/testFlow3d
make
./testFlow3d 21
popd

make -C tests/integration/freeSurfaceTension2d
make -C tests/integration/freeSurfaceTension3d
