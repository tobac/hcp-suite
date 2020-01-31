#!/bin/bash


TASK_LIST="SOCIAL MOTOR LANGUAGE WM GAMBLING RELATIONAL EMOTION"
PARCELLATION_LIST="Glasser Ren RenGlasser"
LIB_DIR="$(dirname $(realpath ${0}))/../lib"
LOG_FILE="/tmp/$(basename ${0}).log"
source "${LIB_DIR}/common.sh"

main() {
  parse_args "${@}"
  parse_common_opts ${ARGS}; log_everything # parse_args sets shifted ARGS for getopts
  ${action}
}

parse_args() {
  # Parse supplied arguments, define basic variables/constants, resolve paths
  # and print usage information

  usage() {
    echo "Usage: ${0} <action> <required arguments for action> [<options>]"
    echo "  Actions:"
    echo "    prepare <task> <IDs file> <parcellation>"
    echo "    design <IDs file> <variables (e.g. \"Age Gender BMI\")> <output file> <CSV table1> [<CSV table2>]"
    echo "    eb <IDs file> <CSV with data> <output file>"
    echo "    run <cope contrast number> <target directory> <IDs file> <task> <parcellation> [-n <number of permutations>] [-o <additional PALM arguments>]"
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
      ARGS="${@}"
    fi
    ;;
  *)
    usage
    ;;
esac

}

## Action functions

prepare() {
  prepare_main() {
    source "$(dirname $(realpath ${0}))/../lib/task_prepare.sh"
    print_summary
    timer start; call_pipeline_script; timer stop Prepare_tfMRI
    timer report
  }

print_summary() {
  echo
  echo "-------------------------------------------------------------"
  echo "Task is ${TASK}."
  echo "Parcellation is ${PARCELLATION}."
  echo "IDs file is ${IDSFILE} -> ${NSUBJECTS} subjects"
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
    echo
  else
    echo
    echo "All subjects prepared successfully."
    echo
  fi
}

prepare_main # Run main action function
}

design() {
  main_design() {
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

create_octave_script() {
  # MATLAB/Octave is incredibly slow for this purpose. picktraits.m takes minutes.
  # R takes ~ 1 sec
  local SCRIPT_OCTAVE="$(create_temp_dir)/build_design.m"
  echo "Script: ${SCRIPT_OCTAVE}"
  cat > ${SCRIPT_OCTAVE} <<EOF
addpath ${LIB_DIR}/HCP/share;
addpath ${FSLDIR}/etc/matlab;
addpath ${HCPPIPEDIR}/global/matlab;

ids_file = "${IDSFILE}";
ids = load(ids_file);
table = "${TABLE}";
outfile = "${OUTFILE_DESIGN}";
EOF
for variable in ${VARIABLES}; do
  echo "printf(\"Picking ${variable}...\n\");" >> ${SCRIPT_OCTAVE}
  echo "${variable} = picktraits(table,{\"${variable}\"},ids,true,'');" >> ${SCRIPT_OCTAVE}
  if [ "${variable}" = "Gender" ]; then # We need to binarise gender
    echo "gender = double(gender == 'M');" >> "${SCRIPT_OCTAVE}"
  fi
done

echo "print(\"Assembling matrix...\n\")" >> ${SCRIPT_OCTAVE}
echo "M = [${VARIABLES} ones(size(${variable}))];" >> ${SCRIPT_OCTAVE}
echo "palm_vestwrite(sprintf('%s',outfile),M);" >> ${SCRIPT_OCTAVE}
}

create_r_script() {
  local TMP_DIR="$(create_temp_dir)"
  SCRIPT_R="${TMP_DIR}"/merge_tables.r
  echo "R script: ${SCRIPT_R}"
  TABLE="${TMP_DIR}"/merged.csv
  cat > "${SCRIPT_R}" <<EOF

data1 <- read.csv("${TABLES[0]}")
data2 <- read.csv("${TABLES[1]}")
cat("Merging tables...\n")
merged <- merge(data1, data2, by="Subject")
write.csv(merged, file = "${TABLE}")

ids <- scan("${IDSFILE}")
vars <- strsplit("${VARIABLES}", split = " ")[[1]] # Convert string into character vector
cat("Picking variables and subjects...\n")
subset <- merged[merged\$Subject %in% ids, ] # First get IDs of interest
subset <- subset[vars] # Then get columns of interest
if ("Gender" %in% names(subset)) {
  subset\$Gender <- ifelse(subset\$Gender == 'M',1,0) # We need to binarise gender
}
subset\$Group <- 1 # Add group column
subset <- sapply(subset, as.numeric) # Convert all columns to numeric

outfile = "${OUTFILE_DESIGN}"
cat("Writing design matrix...\n")
write(paste("/NumWaves", ncol(subset)), outfile)
write(paste("/NumPoints", nrow(subset)), outfile, append = TRUE)
write("", outfile, append = TRUE)
write("/Matrix", outfile, append = TRUE)
#options(scipen = -10) # We want scientific notation in our outfile
write.table(subset, file = outfile, row.names = FALSE, col.names = FALSE, append = TRUE)
#options(scipen = 0)
EOF
}

build_design_with_r() {
  # Keep it simple. Instead of dynamically creating an RScript each time (see above)
  # or feeding arguments to one, we use a static RScript that is very easy to study 
  # and understand and put data into predefined places. This also underlines the role 
  # of this script (and not helper scripts), which is to conveniently parse commonly
  # used arguments/options

  local TMP_DIR="$(create_temp_dir)"
  log cp "${IDSFILE}" "${TMP_DIR}"/ids.file
  log echo "${VARIABLES}" > "${TMP_DIR}"/vars.file
  log cp "${TABLES[0]}" "${TMP_DIR}"/data1.csv
  if [ -n "${TABLES[1]}" ]; then
    log cp "${TABLES[1]}" "${TMP_DIR}"/data2.csv
  fi

  log cd "${TMP_DIR}"
  log Rscript "${LIB_DIR}"/build_design.r
  cd - >/dev/null
  log cp "${TMP_DIR}"/design.mat "${OUTFILE_DESIGN}"
}

main_design
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

main_run # Run main action function

}

# Call main function
main "${@}"
