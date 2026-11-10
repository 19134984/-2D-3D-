#!/bin/bash

PROCESSORS=1
INPUT_FILE=""
NCELLS=1
NRUNS=5

while [[ $# -gt 0 ]]; do
  case $1 in
      --size)
	NCELLS=$2
	shift
	;;
      --runs)
	NRUNS=$2
	shift
	;;
    *)
      INPUT_FILE="$1"
      shift
      ;;
  esac
done


mkdir -p tmp_results

grep -v '^\s*$' "$INPUT_FILE" | while read -r dynamics; do
    ./run_dynamics.sh --runs $NRUNS --size $NCELLS "$dynamics"
done

# Merge all the clean JSON files
jq -s '.' tmp_results/*.json > "${NCELLS}.json"


# Clean up
rm -rf tmp_results

