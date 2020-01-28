#!/bin/bash

source "$(dirname $(realpath ${0}))/../lib/common.sh"
LOG_FILE="/tmp/$(basename ${0}).log"

main() {
  parse_args "${@}"
  parse_common_opts "${ARGS}"; log_everything # parse_args sets shifted ARGS for getopts
  ${action}
}

parse_args() {
  # Parse supplied arguments, define basic variables/constants, resolve paths
  # and print usage information
  TASK_LIST="SOCIAL MOTOR LANGUAGE WM GAMBLING RELATIONAL EMOTION"
  PARCELLATION_LIST="Glasser Ren RenGlasser"

  usage() {
    echo "Usage: ${0} <prepare|run> ..."
    echo "  ${0} prepare <task> <IDs file> <parcellations> [-l <log file>] [-e (do not exit on error)]"
    echo "  ${0} run <cope contrast number> <target directory> <IDs file> <task> <parcellation> [-n <number of permutations>] [-l <log file>] [-e (do not exit on error)] [-o <additional PALM arguments>]"
    exit 1
  }

parse_task() {
  # Check if supplied task is a valid one
  readonly TASK="${1}"
  if ! [[ ${TASK_LIST} =~ (^|[[:space:]])"${TASK}"($|[[:space:]]) ]]; then
    echo "ERROR: Task must be one of the following: ${TASK_LIST}"
    exit 1
  fi
}

parse_parcellations() {
  readonly PARCELLATION="${1}"
  readonly PARCELLATION_FILE="${PARCELLATION_DIR}/${PARCELLATION}.dlabel.nii"
  if ! [[ ${PARCELLATION_LIST} =~ (^|[[:space:]])"${PARCELLATION}"($|[[:space:]]) ]]; then
    echo "ERROR: Parcellation must be one of the following: ${PARCELLATION_LIST}"
    exit 1
  fi
  if ! [ -r "${PARCELLATION_FILE}" ]; then
    echo "ERROR: Parcellation file ${PARCELLATION_FILE} does not exist or is unreadable."
    exit 1
  fi
}

case "${1}" in
  prepare)
    shift
    if [ $# -lt 3 ]; then
      usage
    else
      parse_task "${1}"
      shift # We need to use shift here to get rid of positional arguments for getopts
      parse_idsfile "${1}"
      shift
      parse_parcellations "${1}"
      shift
      action=prepare
      ARGS="${@}"
    fi
    ;;
  run)
    shift
    if [ $# -lt 4 ]; then
      usage
    else
      COPECON="${1}"; shift
      TARGETDIR="$(realpath ${1})"; shift
      TARGET="$(basename ${TARGETDIR})"
      parse_idsfile "${1}"; shift
      parse_task "${1}"; shift
      parse_parcellations "${1}"; shift
      resolve_palm # from common.sh
      action=run
    fi
    ;;
  *)
    usage
    ;;
esac

}

prepare() {
  prepare_main() {
    source "$(dirname $(realpath ${0}))/../lib/task_prepare.sh"
    print_summary
    call_pipeline_script
  }

print_summary() {
  echo
  echo "-------------------------------------------------------------"
  echo "Task is ${TASK}."
  echo "IDs file is ${IDSFILE} -> ${NSUBJECTS} subjects"
  echo "Parcellation is ${PARCELLATION}."
  echo "Log file is ${LOG_FILE}."
  echo "Error handling is ${ERROR_HANDLING}."
  echo "-------------------------------------------------------------"
  echo
  echo -n "Starting in..."
  for n in $(seq 3 -1 1); do echo -n " ${n}"; sleep 1; done
  echo
}

call_pipeline_script() {
  TASKDIR="tfMRI_${TASK}"
  RUNS="${TASKDIR}_RL@${TASKDIR}_LR" # Runs are separated by '@'
  local n=1
  for subject in ${SUBJECTS}; do
    echo -n "Processing subject ${subject} (${n}/${NSUBJECTS})... "
    log ${HCPPIPEDIR}/TaskfMRIAnalysis/TaskfMRIAnalysis.sh \
      --path="${DATADIR}" \
      --subject="${subject}" \
      --lvl1tasks="${RUNS}" \
      --lvl1fsfs="${RUNS}" \
      --lvl2task="${TASKDIR}" \
      --lvl2fsf="${TASKDIR}" \
      --lowresmesh="${LOW_RES_MESH}" \
      --grayordinatesres="${GRAYORDINATES_RESOLUTION}" \
      --origsmoothingFWHM="${ORIGINAL_SMOOTHING_FWHM}" \
      --confound="${CONFOUND}" \
      --finalsmoothingFWHM="${SMOOTHING_FWHM}" \
      --temporalfilter="${TEMPORAL_FILTER}" \
      --vba="${VOLUME_BASED_PROCESSING}" \
      --regname="${REG_NAME}" \
      --parcellation="${PARCELLATION}" \
      --parcellationfile="${PARCELLATION_FILE}" >&3 # Only log to file not to stdout (stderr -> stderr)

    if [ ${?} -eq 0 ]; then # Only makes sense when exit on error is switched off (-e)
      echo "OK"
    else
      echo "FAILED"
      failed_subjects="${failed_subjects} ${subject}"
    fi  
    n=$((n + 1))
  done

  if [ -n "${failed_subjects}" ]; then
    echo
    echo "List of failed subjects:${failed_subjects}" # First character is always space
  fi
}

prepare_main
}

run() {

  main_run() {
    check_already_done
    print_summary
    create_skeleton
    timer start; prepare_for_palm; timer stop PALM_preprocessing
    timer start; run_palm; timer stop PALM
    timer report
  }

print_summary() {
  echo
  echo "-------------------------------------------------------------"
  echo "Task is ${TASK} with COPE number ${COPECON}."
  echo "Target is ${TARGET}."
  echo "Parcellation is ${PARCELLATION}."
  echo "PALM is at ${PALM}; using $NPERM mutations."
  echo "Error handling is ${ERROR_HANDLING}."
  echo "-------------------------------------------------------------"
  echo
  echo -n "Starting in... "
  for n in $(seq 3 -1 1); do echo -n "${n} "; sleep 1; done
  echo
}

prepare_for_palm() {
  # Assemble palm input file
  CIFTIARGS=""
  for subj in ${SUBJECTS}; do
    CIFTIARGS="${CIFTIARGS} -cifti ${DATADIR}/${subj}/MNINonLinear/Results/tfMRI_${TASK}/tfMRI_${TASK}_hp200_s2_level2_${PARCELLATION}.feat/ParcellatedStats/cope${COPECON}.feat/cope1.ptseries.nii -column 1"
  done

  echo -n "Merging input files ..."
  log wb_command -cifti-merge "${TARGETDIR}"/data.ptseries.nii ${CIFTIARGS} && echo OK
}

run_palm() {
  echo "Starting analysis with PALM ..."
  # Running PALM
  log ${PALM} -i "${TARGETDIR}"/data.ptseries.nii -transposedata -d "${TARGETDIR}"/../design/design_${TASK}.mat -t "${TARGETDIR}"/../design/design.con -eb "${TARGETDIR}"/../design/eb_${TASK}.csv -o "${TARGETDIR}"/palm/results -n ${NPERM} -logp -corrcon -accel tail -nouncorrected ${ADDITIONAL_PALM_ARGS}

  for contrast in $(ls -1 "${TARGETDIR}"/palm/results_*cfwep_c*nii | sed -r 's/.*_c([0-9]*)\..*/\1/g'); do 
    log cp "${TARGETDIR}"/palm/results_dat_tstat_cfwep_c${contrast}.pscalar.nii "${TARGETDIR}"/../00-Results/${TARGET}-c${contrast}_${PARCELLATION}.pscalar.nii
  done
}

main_run

}

# Call main function
main "${@}"
