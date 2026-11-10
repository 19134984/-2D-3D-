#!/usr/bin/env bash

current_dir=$(dirname "${BASH_SOURCE[0]}")
source "$current_dir/examples_make_run_check.sh"

# Build core
make

# Array of folders to test
folders=("examples/particles/settlingCube3d")

# Iterate over folders and run the function
for folder in "${folders[@]}"; do
  run_command=$(generate_run_command_parallel "$folder")
  run_make_and_check_output "$run_command" "$folder"
done
