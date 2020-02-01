#!/bin/bash

if [ $# -lt 2 ]; then
  echo
  echo "Usage: ${0} <output cifti> <input cifti 1> <input cifti 2> ... <input cifti n>"
  echo
  exit 1
fi

OUTPUT="${1}"; shift
INPUTS="${@}"

TMPFILE="/tmp/ciftimergeandlabel-${RANDOM}.map"

for input in ${INPUTS}; do
 name=${input//.?scalar.nii/} # Remove suffix
 basename ${name} >> "${TMPFILE}"
done

wb_shortcuts -cifti-concatenate "${OUTPUT}" ${INPUTS}
wb_command -cifti-change-mapping "${OUTPUT}" ROW "${OUTPUT}" -scalar -name-file "${TMPFILE}"
