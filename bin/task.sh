#!/bin/bash

DEPENDENT=true # Tell action functions they are not run standalone
LIB_DIR="$(dirname $(realpath ${0}))/../lib"
source "${LIB_DIR}/common.sh"

main() {
  case ${1} in
    prepare|design|eb|analyse)
      action=${1}
      ;;
    *)
      usage
      ;;
  esac
  shift
  
  parse_args "${@}"
  parse_common_opts ${ARGS} # parse_args sets shifted ARGS for getopts
  log_everything
 
  case ${action} in # Check if all required arguments are there
    prepare)
      if [ ${#mandatory_prepare[@]} -ne 3 ]; then usage; exit 1; fi
      ;;
    design)
      if [ ${#mandatory_design[@]} -ne 4 ]; then usage; exit 1; fi
      ;;
    eb)
      if [ ${#mandatory_eb[@]} -ne 3 ]; then usage; exit 1; fi
      ;;
    analyse)
      if [ ${#mandatory_analyse[@]} -ne 5 ]; then usage; exit 1; fi
      ;;
    *)
      exit # FOR NOW
      ;;
  esac 
  
  source "${LIB_DIR}/task_${action}.sh" # Run action script
}

parse_args() {
  local OPTIND opt OPTARG
  while getopts ":c:d:f:o:p:s:t:V:" opt ${@}; do
    case ${opt} in
      c)
        if [ ${OPTARG} -ne ${OPTARG} ]; then
          echo "Cope contrast must be a number."
          exit 1
        fi
        COPECON=${OPTARG}
        mandatory_analyse+=(1)
        ;;
      d)
        readonly TARGETDIR=$(realpath ${OPTARG})
        readonly TARGET=$(basename ${OPTARG})
        mandatory_analyse+=(1)
        ;;
      f)
        for file in "${OPTARG}"; do
          if ! [ -r "${file}" ]; then
            echo "ERROR: Specified CSV file is not readable."
            exit 1
          fi
          TABLE=$(realpath "${file}")
          TABLES+=($(realpath "${file}"))
        done
        if [ -z ${tables_counted} ]; then
          mandatory_design+=(1)
          mandatory_eb+=(1)
          tables_counted=1
        fi
        ;;     
      o)
        readonly OUTPUT=$(realpath ${OPTARG})
        mandatory_design+=(1)
        mandatory_eb+=(1)
        ;;
      p)
        readonly PARCELLATION="${OPTARG}"
        readonly PARCELLATION_FILE="${PARCELLATION_DIR}/${PARCELLATION}.dlabel.nii"
        if ! [[ ${PARCELLATION_LIST} =~ (^|[[:space:]])"${PARCELLATION}"($|[[:space:]]) ]]; then
          echo "ERROR: Parcellation must be one of the following: ${PARCELLATION_LIST}"
          exit 1
        fi
        if ! [ -r "${PARCELLATION_FILE}" ]; then
          echo "ERROR: Parcellation file ${PARCELLATION_FILE} does not exist or is unreadable."
          exit 1
        fi
        mandatory_prepare+=(1)
        mandatory_analyse+=(1)
        ;;
      s)
        get_subjects "${OPTARG}" # sets IDSFILE, SUBJECTS and NSUBJECTS
        mandatory_prepare+=(1)
        mandatory_design+=(1)
        mandatory_eb+=(1)
        mandatory_analyse+=(1)
        ;;
      t)
        readonly TASK="${OPTARG}"
        if ! [[ ${TASK_LIST} =~ (^|[[:space:]])"${TASK}"($|[[:space:]]) ]]; then
          echo "ERROR: Task must be one of the following: ${TASK_LIST}"
          echo "  Task given: ${TASK}"
          exit 1
        fi
        mandatory_prepare+=(1)
        mandatory_analyse+=(1)
        ;;
      V)
        VARIABLES="${VARIABLES} ${OPTARG}"
        if [ -z ${variables_counted} ]; then mandatory_design+=(1) && variables_counted=1; fi
        ;;
      :)
        echo "ERROR: Invalid option: -${opt} requires an argument." 1>&2
        exit 1
        ;;
    esac
  done
  [ ${OPTIND} -gt 1 ] && shift $((OPTIND - 2)) # Prevent error when < 2 arguments are passed
  ARGS="${@}"
}

usage() {
  echo "Usage: ${0} <action> [action-specific arguments] [common options]"
  echo "  Actions:"
  echo "    prepare"
  echo "      -s <arg>: specify file with subject IDs or a space-separated list of subjects"
  echo "      -t <arg>: specify a task"
  echo "      -p <arg>: specify a parcellation"
  echo "      You might want to consider using -e to prevent exiting when one subject fails."
  echo "    design"
  echo "      -s <arg>: specify file with subject IDs or a space-separated list of subjects"
  echo "      -o <arg>: specify an output file (e.g. design.mat)"
  echo "      -V <arg>: specify variables (e.g. \"Age_in_Yrs Gender BMI\"; repeatable)"
  echo "      -f <arg>: specify CSV file containing the specified variables (repeatable; e.g."
  echo "                \"-f restricted.csv -f unrestricted.csv\")"
  echo "    eb"
  echo "      -s <arg>: specify file with subject IDs or a space-separated list of subjects"
  echo "      -o <arg>: specify an output file (e.g. eb.csv)"
  echo "      -f <arg>: specify CSV file containing the specified variables (repeatable; e.g."
  echo "                \"-f restricted.csv -f unrestricted.csv\")"
  echo "    analyse"
  echo "      -s <arg>: specify file with subject IDs or a space-separated list of subjects"
  echo "      -t <arg>: specify a task"
  echo "      -c <arg>: specify a COPE number"
  echo "      -p <arg>: specify a parcellation"
  echo "      -d <arg>: specify target directory"
  echo "      [-n <arg>]: number of permutations (default: 5000)"
  echo "      [-o <arg>]: additional PALM arguments"
  echo "      Specify design, contrast and EB via additional PALM arguments or place them in"
  echo "      <target directory>/../ as design_<task>.mat, design.con and eb_<task>.csv"
  echo
  echo "  Common options:"
  echo '    [-l <log file>]: specify log file instead of default /tmp/${0}-${RANDOM}.log'
  echo "    [-e]: turn off error handling (do not exit on error)"
  exit 1
}

# Call main function
main "${@}"
