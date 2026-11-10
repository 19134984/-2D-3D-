#!/usr/bin/env bash

cd external/voreen/frontend
nix --extra-experimental-features nix-command --extra-experimental-features flakes build .#voreen
