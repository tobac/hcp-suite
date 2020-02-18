#!/usr/bin/python3
import sys
import numpy as np
import nibabel as nib
from sklearn.preprocessing import StandardScaler
import math
from roi_connectome import *
import subprocess
import tempfile
import os

def get_nimg_data(fname):
    data = nib.load(fname).get_fdata()
    return data

def cifti_convert_to_nifti(input_fname, output_fname):
    subprocess.run(["wb_command", "-cifti-convert", "-to-nifti", input_fname, output_fname], capture_output=True)
    
def run_film_gls(input_fname, output_dir, design_fname, contrast_fname):
    subprocess.run(["film_gls", "--in={}".format(input_fname), "--rn={}".format(os.path.join(output_dir, "film_gls_output")),
                    "--pd={}".format(design_fname), "--con={}".format(contrast_fname), "--thr=1", "--mode=volumetric", "-v"])
    return os.path.join(output_dir, "film_gls_output", "prewhitened_data.nii.gz")

def prewhiten_cifti(cifti_fname, design_fname, contrast_fname, tmpdir):
    fakenifti_fname = os.path.join(tmpdir.name, "fakenifti.nii.gz")
    # We could use singlesubject_cifti_to_fake_nifti(), but we have to save a temporary file anyway for film_gls and wb_command is certainly more robust
    cifti_convert_to_nifti(cifti_fname, fakenifti_fname)
    prewhitened_data_fname = run_film_gls(fakenifti_fname, tmpdir.name, design_fname, contrast_fname)
    return prewhitened_data_fname

def demean_and_standardize(array):
    scaler = StandardScaler()
    scaler.fit(array)
    demeaned_and_standardized_array = scaler.transform(array)
    return demeaned_and_standardized_array

def get_ev_timeseries(ids, task, runs, parcellation, ev_files, data_dir='/home/tobac/HCP/S1200/', tr=0.72):
    ev_data_dict = {}
    n = 1
    for id in ids:
        ev_data_dict[id] = [] # Create subject-specific list to store EV data in
        T = 0
        print("Processing {}... ({}/{} | {} %)".format(id, n, len(ids), round(n/len(ids)*100, 2)))
        for run in runs:
            tmpdir = tempfile.TemporaryDirectory()
            run_path = os.path.join(data_dir, str(id), "MNINonLinear", "Results", "tfMRI_{}_{}".format(task, run))
            cifti_fname = os.path.join(run_path, "tfMRI_{}_{}_Atlas_hp200_s2_{}.ptseries.nii".format(task, run, parcellation))
            feat_dir = os.path.join(run_path, "tfMRI_{}_{}_hp200_s2_level1_{}.feat".format(task, run, parcellation))
            design_fname = os.path.join(feat_dir, "design.mat")
            contrast_fname = os.path.join(feat_dir, "design.con")
            
            # Preprocessing: Prewhiten
            prewhitened_data_fname = prewhiten_cifti(cifti_fname, design_fname, contrast_fname, tmpdir) 
            nifti_data = get_nimg_data(prewhitened_data_fname)
            
            
            # Downstream operations will be easier if we use timepoints:parcels shape instead of parcels:1:1:timepoints
            subject_data = nifti_data[:, 0, 0, :] # Fancy indexing gets rid of the superfluous middle dimensions -> parcels:timepoints
            subject_data = np.swapaxes(subject_data, 0, 1) # Finally swap axes -> timepoints:parcels
            
            # Demean and standardize (aka "fit-transform") the data to deal with multiple runs
            subject_data_std = demean_and_standardize(subject_data)
            
            for ev_file in ev_files:
                ev_fname = os.path.join(run_path, "EVs", ev_file)
                evs = np.loadtxt(ev_fname)
                e = 1
                for ev in evs:
                    start = math.ceil((ev[0]/tr)+1) # We can afford to lose one timepoint at the beginning, it might even improve SNR
                    end = math.floor((((ev[0]+ev[1])/tr)+1)) # We can afford to lose one timepoint at the end, it might even improve SNR
                    ev_data = subject_data_std[start:end]
                    print("\tTimepoints in EV part {}/{} of {} (run {}): {}".format(e, len(evs), ev_file, run, len(ev_data)))
                    T+=len(ev_data)
                    for ev_data_item in ev_data: # Don't create nested lists, just append each item to the primary list
                        ev_data_dict[id].append(ev_data_item)
                    e+=1
            tmpdir.cleanup() # Delete temporary directory
        print("\t-> Total timepoints for subject {}: {}".format(id, T))
        n+=1
    return ev_data_dict

def save_data_dict(data_dict):
    n = 1
    N = len(data_dict)
    fname_list = []
    for id, data in data_dict.items():
        fname = "{}.txt".format(id)
        print("Saving text file {}/{} ({} %)".format(n, N, round(n/N*100, 2)), end="\r")
        np.savetxt(fname, data)
        fname_list.append(fname)
        n+=1
        sys.stdout.flush # Needed to update line in console
    return fname_list

def singlesubject_cifti_to_fake_nifti(cifti_data):
    newarray = np.array(np.zeros((cifti_data.shape[1], 1, 1, cifti_data.shape[0])))
    for idx in range(0, len(cifti_data)): 
        for idy in range(0, len(cifti_data[idx])): 
            newarray[idy][0][0][idx] = cifti_data[idx][idy]
    affine = [[1, 0, 0, 0],[0, 1, 0, 0],[0, 0, 1, 0],[0, 0, 0, 1]]
    img = nib.nifti1.Nifti1Image(newarray, affine=affine)
    return img