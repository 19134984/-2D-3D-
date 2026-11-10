#!/usr/bin/env bash

tar xfv olb-master.tgz

set -eou pipefail

pushd olb-master
make
make samples
popd
