#!/bin/bash

main() {
  usage $@
  extract_ic
  separate_cifti
  get_clusters
}

extract_ic() {
  if [ $(wb_command -file-information "${ICFILE}" -only-number-of-maps) -gt 1 ]; then
    echo -n "Extracting IC #${IC} from file ${ICFILE} ..."
    log wb_shortcuts -cifti-concatenate -map ${IC} ${TMPDIR}/ic${IC}.dscalar.nii "${ICFILE}"
    echo " OK"
  else
    echo "Assuming that IC in question is already isolated in file ${ICFILE}."
    cp "${ICFILE}" ${TMPDIR}/ic${IC}.dscalar.nii
  fi
}

separate_cifti() {
  echo -n "Extracting left and right cerebellum ..."
  for side in left right; do 
    # First extract each cerebellar side as volume (NIFTI)
    log wb_command -cifti-separate ${TMPDIR}/ic${IC}.dscalar.nii COLUMN -volume CEREBELLUM_$(echo ${side} | awk '{print toupper($0)}') ${TMPDIR}/ic${IC}-cerebellum-${side}.nii
    # Then recreate CIFTI file with volume as only content
    log wb_command -cifti-create-dense-from-template ${TMPDIR}/ic${IC}.dscalar.nii ${TMPDIR}/ic${IC}-cerebellum-${side}.dscalar.nii -volume CEREBELLUM_$(echo ${side} | awk '{print toupper($0)}') ${TMPDIR}/ic${IC}-cerebellum-${side}.nii
    set +v
  done
  echo " OK"
}

get_clusters() {
  for side in left right; do
    # 1. Actual cluster finding
    echo -n "Finding clusters in ${side} cerebellum ..."
    log wb_command -cifti-find-clusters ${TMPDIR}/ic${IC}-cerebellum-${side}.dscalar.nii 0 0 ${VTHR} ${STHR} COLUMN ${TMPDIR}/ic${IC}-clusters-cerebellum-${side}-clusters-positive.dscalar.nii -left-surface ${SRCDIR}/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii -right-surface S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii
    log wb_command -cifti-find-clusters ${TMPDIR}/ic${IC}-cerebellum-${side}.dscalar.nii 0 0 -${VTHR} ${STHR} COLUMN ${TMPDIR}/ic${IC}-clusters-cerebellum-${side}-clusters-negative.dscalar.nii -left-surface ${SRCDIR}/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii -right-surface S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii -less-than
    echo " OK"
    # 2. Convert clusters' dense scalar to dense label
    for direction in positive negative; do
      log wb_command -cifti-label-import ${TMPDIR}/ic${IC}-clusters-cerebellum-${side}-clusters-${direction}.dscalar.nii '' ${TMPDIR}/ic${IC}-clusters-cerebellum-${side}-clusters-${direction}.dlabel.nii
      n_clusters=$(($(wb_command -file-information ${TMPDIR}/ic${IC}-clusters-cerebellum-${side}-clusters-${direction}.dlabel.nii -only-cifti-xml | grep -c "Label Key")-1))
      echo "  Number of ${direction} clusters in ${side} cerebellum: ${n_clusters}"
      # 3. Convert each label to a ROI file
      for cluster in $(seq -w 01 ${n_clusters}); do
        basename="ic${IC}-cerebellum_${side}-${VTHR}_${STHR}-${direction}_cluster_${cluster}"
        test ${n_clusters} -gt 0 && log wb_command -cifti-label-to-roi ${TMPDIR}/ic${IC}-clusters-cerebellum-${side}-clusters-${direction}.dlabel.nii ${TMPDIR}/results/${basename}.dscalar.nii -key $(echo ${cluster} | sed 's/^0//g')
        print_meta
      done
    done
  done
  echo 
  echo 
  echo "Resulting ROIs reside in:"
  echo "  ${TMPDIR}/results"
  echo
}

print_meta() {
  metafile="${TMPDIR}/results/${basename}.meta"
  echo "IC number: ${IC}" > "${metafile}"
  echo "Value threshold: ${VTHR}" >> "${metafile}"
  echo "Size threshold in mm^3: ${STHR}" >> "${metafile}"
  echo >> "${metafile}"
  echo "Relevant commands used:" >> "${metafile}"
  cat "${TMPDIR}"/log >> "${metafile}"
}

usage() {
  if [ $# -lt 3 ]; then
    echo
    echo "Usage: $(basename ${0}) <IC number> <value threshold, e.g. 2*stddev> <size threshold [mm^3]> [IC source file]"
    echo
    echo "$(basename ${0}) assumes that group average source files are located in the same directory as your IC source file."
    exit 1
  fi
  IC=$(printf "%02d\n" ${1})
  VTHR=${2}
  STHR=${3}
  if [ -z ${4} ]; then ICFILE="./melodic_IC.dscalar.nii"; else ICFILE=${4}; fi
  SRCDIR=$(dirname $(realpath ${ICFILE}))
  TMPDIR=/tmp/ic2crbllmroi-${RANDOM}
  mkdir -p ${TMPDIR}/results
}

log() {
  command=("$@")
  echo "${command[@]}" >> "${TMPDIR}"/log
  "${command[@]}"
}

set -e
main $@
