#!/bin/bash

PROCESSORS=6
INPUT_FILE=""
PRECISION=float

while [[ $# -gt 0 ]]; do
  case $1 in
      --precision)
	PRECISION=$2
	shift
	;;
    *)
      INPUT_FILE="$1"
      shift
      ;;
  esac
  shift
done

echo "$PRECISION"

cat $INPUT_FILE | xargs -I {} -P $PROCESSORS ./compile_dynamics.sh --precision $PRECISION "{}"
