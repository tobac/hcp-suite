#!/bin/bash

run() {
  print_summary
  create_octave_script
  octave -qf "${SCRIPT_OCTAVE}"

}

print_summary() {
  echo
  echo "-------------------------------------------------------------"
  echo "IDs file is ${IDSFILE} -> ${NSUBJECTS} subjects"
  echo "Table is ${TABLE}."
  echo "Output file is ${OUTPUT}."
  echo "Log file is ${LOG_FILE}."
  echo "Error handling is ${ERROR_HANDLING}."
  echo "-------------------------------------------------------------"
  echo
  echo -n "Starting in..."
  for n in $(seq 3 -1 1); do echo -n " ${n}"; sleep 1; done
  echo
}

create_octave_script() {
  local TMP_DIR="$(create_temp_dir)"
  SCRIPT_OCTAVE="${TMP_DIR}"/build_eb.m
  cat > "${SCRIPT_OCTAVE}" <<EOF
addpath ${LIB_DIR}/HCP/share;

ids = load("${IDSFILE}");
printf("Building exchangeability blocks...\n")
EB = hcp2blocks("${TABLE}","${OUTPUT}",true,ids);
EOF
}

run # Call main function
