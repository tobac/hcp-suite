#!/bin/bash

DEPENDENT=true # Tell action functions they are not run standalone
LIB_DIR="$(dirname $(realpath ${0}))/../lib"
source "${LIB_DIR}/common.sh"

main() {
  parse_args "${@}"
  parse_common_opts ${ARGS} # parse_args sets shifted ARGS for getopts
  log_everything
  source "${LIB_DIR}/task_${action}.sh" # Run action script
}

parse_args() {
  # Parse supplied arguments, define basic variables/constants, resolve paths
  # and print usage information

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
    design)
      shift
      if [ $# -lt 4 ]; then
        usage
      else
        parse_idsfile "${1}"
        shift
        VARIABLES="${1}"
        shift
        OUTFILE_DESIGN="${1}"
        shift
        while ! [[ "${1}" =~ ^- ]] && ! [ -z "${1}" ]; do
          TABLES+=("${1}") # Add element to array
          shift
        done
        action=design
        ARGS="${@}"
      fi
      ;;
    eb)
      shift
      if [ $# -lt 3 ]; then
        usage
      else
        parse_idsfile "${1}"
        shift
        TABLE="${1}"
        shift
        OUTFILE_EB="${1}"
      fi
      action="source ${LIB_DIR}/task_eb.sh"
      ARGS="${@}"
      ;;
    analyse) 
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
        action=analyse
        ARGS="${@}"
      fi
      ;;
    *)
      usage
      ;;
  esac
}

usage() {
  echo "Usage: ${0} <action> <required arguments for action> [<options>]"
  echo "  Actions:"
  echo "    prepare <task> <IDs file> <parcellation>"
  echo "    design <IDs file> <variables (e.g. \"Age Gender BMI\")> <output file> <CSV table1> [<CSV table2>]"
  echo "    eb <IDs file> <CSV with data> <output file>"
  echo "    analyse <cope contrast number> <target directory> <IDs file> <task> <parcellation> [-n <number of permutations>] [-o <additional PALM arguments>]"
  echo "      Specify design, contrast and EB via additional PALM arguments or place them in"
  echo "      <target directory>/../ as design_<task>.mat, design.con and eb_<task>.csv"
  echo
  echo "  Common options:"
  echo '    [-l <log file>]: specify log file instead of default /tmp/${0}-${RANDOM}.log'
  echo "    [-e]: turn off error handling (do not exit on error)"
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

# Call main function
main "${@}"
