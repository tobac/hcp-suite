#!/bin/bash 
# This is a combination of HCP Suite specific options (adopted from
# TaskfMRIAnalysisBatch.sh) and SetUpHCPPipeline.sh from 
# ${HCPPIPEDIR}/Examples/Scripts

## adopted from TaskfMRIAnalysisBatch.sh

export CARET7DIR="/usr/bin" # Path to wb_command
SMOOTHING_FWHM="2" # Smoothing in mm. 2 = ORIGINAL_SMOOTHING_FWHM -> no more smoothing than already added by HCP minimal preprocessing pipelines
ORIGINAL_SMOOTHING_FWHM="2" # 2 mm if using HCP minimal preprocessing pipeline outputes
LOW_RES_MESH="32" # 32 if using HCP minimal preprocessing pipeline outputs
GRAYORDINATES_RESOLUTION="2" # 2 mm if using HCP minimal preprocessing pipeline outputs
CONFOUND="NONE" # File located in ${SubjectID}/MNINonLinear/Results/${fMRIName} or NONE
TEMPORAL_FILTER="200" # Use 2000 for linear detrend, 200 is default for HCP task fMRI
VOLUME_BASED_PROCESSING="NO" # YES or NO. CAUTION: Only use YES if you want unconstrained volumetric blurring of your data, otherwise set to NO for faster, less biased, and more senstive processing (grayordinates results do not use unconstrained volumetric blurring and are always produced).
REG_NAME="NONE" # Use NONE to use the default surface registration


if [ -z "${FSLDIR}" ]; then
  echo '${FSLDIR} is not set but required.' 1>&2
  exit 1
fi

## from SetUPHCPPipeline.sh

# Let FreeSurfer know what version of FSL to use
# FreeSurfer uses FSL_DIR instead of FSLDIR to determine the FSL version
export FSL_DIR="${FSLDIR}"

# Set up FreeSurfer (if not already done so in the running environment)
# Uncomment the following 2 lines (remove the leading #) and correct the FREESURFER_HOME setting for your setup
#export FREESURFER_HOME=/usr/local/bin/freesurfer
#source ${FREESURFER_HOME}/SetUpFreeSurfer.sh > /dev/null 2>&1

# Set up specific environment variables for the HCP Pipeline
export HCPPIPEDIR="${HCP}"/HCPpipelines

export HCPPIPEDIR_Templates=${HCPPIPEDIR}/global/templates
export HCPPIPEDIR_Bin=${HCPPIPEDIR}/global/binaries
export HCPPIPEDIR_Config=${HCPPIPEDIR}/global/config

export HCPPIPEDIR_PreFS=${HCPPIPEDIR}/PreFreeSurfer/scripts
export HCPPIPEDIR_FS=${HCPPIPEDIR}/FreeSurfer/scripts
export HCPPIPEDIR_PostFS=${HCPPIPEDIR}/PostFreeSurfer/scripts
export HCPPIPEDIR_fMRISurf=${HCPPIPEDIR}/fMRISurface/scripts
export HCPPIPEDIR_fMRIVol=${HCPPIPEDIR}/fMRIVolume/scripts
export HCPPIPEDIR_tfMRI=${HCPPIPEDIR}/tfMRI/scripts
export HCPPIPEDIR_dMRI=${HCPPIPEDIR}/DiffusionPreprocessing/scripts
export HCPPIPEDIR_dMRITract=${HCPPIPEDIR}/DiffusionTractography/scripts
export HCPPIPEDIR_Global=${HCPPIPEDIR}/global/scripts
export HCPPIPEDIR_tfMRIAnalysis=${HCPPIPEDIR}/TaskfMRIAnalysis/scripts

#try to reduce strangeness from locale and other environment settings
export LC_ALL=C
export LANGUAGE=C
#POSIXLY_CORRECT currently gets set by many versions of fsl_sub, unfortunately, but at least don't pass it in if the user has it set in their usual environment
unset POSIXLY_CORRECT
