#!/bin/bash

# For standalone use, the following variables need to be defined:
#IDSFILE="subjects.txt"
#TABLES=("restricted.csv" "unrestricted.csv") # As array
#OUTFILE_DESIGN="design.mat"
#VARIABLES="Age_in_Yrs Gender"

main() {
  print_summary
  if [ ${#TABLES[@]} -eq 1 ] || [ ${#TABLES[@]} -eq 2 ]; then                                
    build_design_with_r
  else
    echo "Specify either one or two CSV files."                                              
    exit 1
  fi
}

print_summary() {
  echo
  echo "-------------------------------------------------------------"                         
  echo "IDs file is ${IDSFILE}"
  echo "Data CSV is/are ${TABLES[@]}."
  echo "Output file is ${OUTFILE_DESIGN}."
  echo "Variables are:"
  echo "  ${VARIABLES}"
  echo "Log file is ${LOG_FILE}."
  echo "Error handling is ${ERROR_HANDLING}."
  echo "-------------------------------------------------------------"                         
  echo
  echo -n "Starting in... "
  for n in $(seq 3 -1 1); do echo -n "${n} "; sleep 1; done                                    
  echo
}

build_design_with_r() {
  # Keep it simple. Instead of dynamically creating an RScript each time
  # or feeding arguments to one, we use a static RScript that is very easy to study 
  # and understand and put data into predefined places. This also underlines the role 
  # of the wrapper script (and not helper scripts), which is to conveniently parse 
  # commonly used arguments/options

  local TMP_DIR="$(create_temp_dir)"
  cp "${IDSFILE}" "${TMP_DIR}"/ids.file
  echo "${VARIABLES}" > "${TMP_DIR}"/vars.file
  cp "${TABLES[0]}" "${TMP_DIR}"/data1.csv
  if [ -n "${TABLES[1]}" ]; then
    cp "${TABLES[1]}" "${TMP_DIR}"/data2.csv
  fi

  cd "${TMP_DIR}"
  Rscript "${LIB_DIR}"/build_design.r
  cd - >/dev/null
  cp "${TMP_DIR}"/design.mat "${OUTFILE_DESIGN}"
}

main
