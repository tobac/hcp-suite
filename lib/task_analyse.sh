#!/bin/bash

main() {
  check_already_done
  resolve_palm
  print_summary
  create_skeleton
  timer start; prepare_for_palm; timer stop PALM_preprocessing
  timer start; run_palm; timer stop PALM
  timer report
}

print_summary() {
  echo
  echo "------------------------------------------------------------------"
  echo "Task is ${TASK} with COPE number ${COPECON}."
  echo "Target is ${TARGET}."
  echo "Parcellation is ${PARCELLATION}."
  echo "PALM is at ${PALM}; using $NPERM mutations."
  echo "Log file is ${LOG_FILE}."
  echo "Error handling is ${ERROR_HANDLING}."
  echo "------------------------------------------------------------------"
  echo
  echo -n "Starting in... "
  for n in $(seq 3 -1 1); do echo -n "${n} "; sleep 1; done
  echo
}

prepare_for_palm() {
  # Assemble palm input file
  local PATH_TO_COPE_CIFTI="MNINonLinear/Results/tfMRI_${TASK}/tfMRI_${TASK}_hp200_s2_level2_${PARCELLATION}.feat/ParcellatedStats/cope${COPECON}.feat/cope1.ptseries.nii"
  local ciftiargs=() # Create a local array to hold wb_command arguments in
  for subject in ${SUBJECTS}; do
    ciftiargs+="-cifti ${DATADIR}/${subject}/${PATH_TO_COPE_CIFTI} -column 1 " # Add each subject's file argument to array
  done

  echo -n "Merging input files ..."
  log wb_command -cifti-merge "${TARGETDIR}"/data.ptseries.nii ${ciftiargs[@]} && echo OK
}

run_palm() {
  echo "Starting analysis with PALM ..."
  # Running PALM
  log ${PALM} -i "${TARGETDIR}"/data.ptseries.nii -transposedata -d "${TARGETDIR}"/../design/design_${TASK}.mat -t "${TARGETDIR}"/../design/design.con -eb "${TARGETDIR}"/../design/eb_${TASK}.csv -o "${TARGETDIR}"/palm/results -n ${NPERM} -logp -corrcon -accel tail -nouncorrected ${ADDITIONAL_PALM_ARGS}

  for contrast in $(ls -1 "${TARGETDIR}"/palm/results_*cfwep_c*nii | sed -r 's/.*_c([0-9]*)\..*/\1/g'); do
    log cp "${TARGETDIR}"/palm/results_dat_tstat_cfwep_c${contrast}.pscalar.nii "${TARGETDIR}"/../00-Results/${TARGET}-c${contrast}_${PARCELLATION}.pscalar.nii
  done
}

if [ -z ${DEPENDENT} ]; then # Source common.sh if run standalone
  source "$(dirname ${0})"/common.sh
fi

main
