#!/bin/bash

PROCESSORS=2
PRECISION=float

while [[ $# -gt 0 ]]; do
  case $1 in
    --precision)
	PRECISION=$2
	shift
	;;
    *)
      INPUT_FILE=$1
      shift
      ;;
  esac
done


# Introspection disabling automatic CSE 
UNOPT_JSON=$(cat $INPUT_FILE | xargs -I {} -P $PROCESSORS ./inspect_dynamics.sh --precision $PRECISION --disable-automatic-cse "{}" | jq -s '.' )

# Introspection with automatic CSE
OPT_JSON=$(cat $INPUT_FILE | xargs -I {} -P $PROCESSORS ./inspect_dynamics.sh --precision $PRECISION "{}" | jq -s '.' )

jq -s 'add' <(echo "$UNOPT_JSON") <(echo "$OPT_JSON") > dynamics.json
