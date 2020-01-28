#!/bin/bash
# Functions and options used commonly by HCP Suite scripts

## Options
HCP="/home/tobac/HCP"
DATADIR="${HCP}/S1200"
PARCELLATION_DIR="/home/tobac/Projekte/diss/parcellations"
LOG_FILE="/tmp/$(basename ${0})-${RANDOM}.log"
ulimit -n 4096

## Functions
set_error_handling() { # On per default -> see last line in this file
  case ${1} in
    on)
      set -o errtrace
      set -o pipefail
      trap 'exit_code=${?}; last_command=${BASH_COMMAND}; function_name=${FUNCNAME[0]}; onerr' ERR
      ERROR_HANDLING="on"
      ;;
    off)
      trap - ERR # Unsetting the trap makes sense in some cases -> don't exit on error
      ERROR_HANDLING="off"
      ;;
  esac
}

onerr() {
  echo
  echo "ERROR:"
  echo "  Command: ${last_command}"
  echo "  finished with exit code: ${exit_code}" 
  echo "  in function: ${function_name}."
  if [ -n "${LOCKFILE}" ]; then
    [ -f "${LOCKFILE}" ] && rm "${LOCKFILE}"
  fi
  exit 1
}

log() {
  command=("${@}")
  echo "${command[@]}" >&3
  "${command[@]}"
}

parse_common_opts() {
  while getopts ":en::l::o:" opt ${@}; do
    case ${opt} in
      e)
        # Don't exit on error
        set_error_handling off
        ;;
      n)
        # Sets number of permutations for PALM
        NPERM=${OPTARG}
        if ! [[ ${NPERM} =~ ^[0-9]+$ ]]; then
          echo "Number of permutations needs to be an integer." 1>&2
          exit 1
        fi
        ;;
      l)
        LOG_FILE="${OPTARG}"
        touch "${LOG_FILE}" || (echo "ERROR: Unable to create log file ${LOG_FILE}" >&2; exit 1)
        ;;
      o)
        ADDITIONAL_PALM_ARGS="${OPTARGS}"
        ;;
      \?)
        echo "ERROR: Invalid option: -${opt}". 1>&2
        exit 1
        ;;
      :)
        echo "ERROR: Invalid option: -${opt} requires an argument." 1>&2
        exit 1
        ;;
      *)
        echo "No valid option specified."
        ;;
    esac
  done
  shift $((OPTIND - 1))
}

parse_idsfile() {
  readonly IDSFILE=$(realpath "${1}")
  if ! [ -r "${IDSFILE}" ]; then
    echo "IDs file does not exist or is unreadable."
    exit 1
  fi
  SUBJECTS=$(cat "${IDSFILE}")
  NSUBJECTS=$(cat "${IDSFILE}" | wc -l) # Pipe to avoid filename output
}

log_everything() {
  #set -x
  exec > >(tee -a ${LOG_FILE} )
  exec 2> >(tee -a ${LOG_FILE} >&2)
  exec 3> >(cat >> ${LOG_FILE}) # only write to log file, not to stdout
}

timer() {
  # Usage: timer <start|stop|report> <timer_name>
  case $1 in
    start)
      starttime=$(date '+%s')
      ;;
    stop)
      stoptime=$(date '+%s')
      elapsed=$((${stoptime}-${starttime}))
      timer_all="${timer_all} ${2}:${elapsed}"
      echo
      echo "Time elapsed for ${2}: ${elapsed} s"
      echo
      ;;
    report)
      echo
      echo "-------------------------------------------------------------------------------"
      echo "Timer report of $(basename ${0}):"
      for timer in ${timer_all}; do
        echo "${timer}" | awk -F ':' '{print "  " $1 ": " $2 " s"}'
      done
      echo "-------------------------------------------------------------------------------"
      echo
  esac
}

check_already_done() {
  if [ -d "${TARGETDIR}" ]; then
    echo ""
    echo "[x] Target ${TARGET} already analysed, exiting."
    echo ""
    exit 1
  fi
}

create_skeleton() {
  for dir in analysis design palm midthick results; do
    test -d "${TARGETDIR}"/${dir} || mkdir -p "${TARGETDIR}"/${dir}
  done
  touch "${TARGETDIR}"/.lock
  if ! [ -d "${TARGETDIR}/../00-Results" ]; then
    mkdir "${TARGETDIR}/../00-Results"
    cp "${HCP}"/S1200_GroupAvg_v1/S1200.{L,R}.midthickness_MSMAll.32k_fs_LR.surf.gii "${TARGETDIR}"/../00-Results/
    cp "${HCP}"/S1200_GroupAvg_v1/S1200_AverageT2w_restore.nii.gz "${TARGETDIR}"/../00-Results/
  fi
  cp "${0}" "${TARGETDIR}"/
  echo "${0} ${ARGS}" > "${TARGETDIR}"/command.log
  if [ -d design ]; then cp design/* "${TARGETDIR}"/design/; fi
  cp "${IDFILE}" "${TARGETDIR}"
}

create_temp_dir() {
  tmpdir="/tmp/hcpsuite-${RANDOM}"
  mkdir "${tmpdir}"
  return "${tmpdir}"
}

resolve_palm() {
  PALM="${HOME}/Downloads/PALM/palm /usr/bin/palm"
  for palm in ${PALM}; do
    if [ -x "${palm}" ]; then
      PALM="${palm}"
    fi
  done
  if ! [ -x "${PALM}" ]; then
    echo "PALM not found"
    exit 1
  fi
}

set_error_handling on
