#!/usr/bin/env bash

# Function to make, run, and check output of examples
run_make_and_check_output() {
  local run_command="$1"
  local folder="$2"
  local name=$(basename "$folder")

  # Read expected output from files
  required_fixed_strings_file="tests/$folder/expected_output_fixed.txt"
  if [ ! -f "$required_fixed_strings_file" ]; then
    echo "Error: $required_fixed_strings_file not found"
    exit 1
  fi
  required_regexp_file="tests/$folder/expected_output_regexp.txt"
  if [ ! -f "$required_regexp_file" ]; then
    echo "Error: $required_regexp_file not found"
    exit 1
  fi

  # Build and fail on error
  make -C "$folder" || exit 1

  # Create a named pipe
  fifo=$(mktemp -u)
  mkfifo "$fifo"

  # Run the command in the background
  pushd $folder
  $run_command 2>&1 | tee "$fifo" &
  popd

  # Read from the named pipe into a variable
  output=$(cat "$fifo")

  # Remove the named pipe
  rm "$fifo"

  # Check if each required string is present in the output
  error_found=0
  while IFS= read -r str; do
    if ! grep -qF "$str" <<< "$output"; then
      echo "Required string was not found in output for $folder: \"$str\""
      error_found=1
    fi
  done < "$required_fixed_strings_file"
  while IFS= read -r str; do
    if ! grep -qE "$str" <<< "$output"; then
      echo "Required regexp was not found in output for $folder: \"$str\""
      error_found=1
    fi
  done < "$required_regexp_file"

  if [ "$error_found" -eq 1 ]; then
    exit 1
  fi

  echo "The output matches the expectations for $folder."
}

generate_run_command_parallel() {
  local folder="$1"
  local name=$(basename "$folder")

  number_processes=$(($(grep -c ^processor /proc/cpuinfo) / 2))

  local run_command="mpirun -np $number_processes $name"

  echo "$run_command"
}

generate_run_command_sequential() {
  local folder="$1"
  local name=$(basename "$folder")
  local run_command="./$name"
  echo "$run_command"
}
