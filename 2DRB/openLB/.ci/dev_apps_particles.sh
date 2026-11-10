#!/usr/bin/env bash

make && \
make -C apps/dev_Particles/ makeIntegrationtests && \
export MAKEFLAGS="-j1" && \
make -C apps/dev_Particles/ runIntegrationtests
