#!/usr/bin/python3

import sys
import os
import tempfile
import subprocess
import math
import time
import numpy as np
import nibabel as nib
from tempfile import mkdtemp
from nilearn.connectome import ConnectivityMeasure
from sklearn.preprocessing import StandardScaler
from nilearn import plotting
from matplotlib import pyplot as plt
from pandas import read_csv


global timers
timers = [] 
def timer(command, name=""):
  msg = ""
  if command == "tic":
    global tic
    tic = time.perf_counter()
  elif command == "toc":
    toc = time.perf_counter()
    msg = "Time elapsed for {}: {} s".format(name, round(toc-tic, 2))
    print("{}\n".format(msg))
    timers.append(msg)
  elif command == "report":
    print("\n------------------------------------------------------------")
    print("Timer report of {}:".format(os.path.basename(sys.argv[0])))
    for timer in timers:
      print("  {}".format(timer))
    print("------------------------------------------------------------\n")


def load_list_from_file(file_name):
  list = [line.rstrip('\n') for line in open(file_name)]
  return list


def load_time_series(file_list_fname):
  """Load time series from pre-specified file"""
  time_series = []
  file_list = load_list_from_file(file_list_fname)
  n = 1
  n_series = len(file_list)

  timer("tic")
  for file in file_list:
    print("Loading time series from file {} ({}/{} | {} %)".format(file, n, n_series, round((n/n_series)*100, 2)))
    # Using Pandas' read_csv() speeds things up enormously in comparison to np.loadtxt
    ts = read_csv(file, delim_whitespace=True)
    tsnp = np.array(ts) # Convert to numpy array
    time_series.append(tsnp)
    n += 1
  timer("toc", name="loading all time series")

  return time_series, file_list


def save_whole_sample_nifti(matrices_array, output_path, output_fname=None, coi_indices=None, coi_fname=None, clean=True):
  """
  Take an array of matrices (e.g. output from ConnecitivityMeasure()),
  rearrange it to build a multi-subject NIfTI, e.g. for 361 ROIs and 
  999 subjects: (999, 361, 361) -> (361, 361, 1, 999). Default filename is
  output_path/corr_matrices.full.nii.gz and can be overriden with output_file.
  
  If coi_indices are supplied, also save a COI matrix at 
  output_path/coi_matrices.nii.gz (default) or coi_fname.
  
  NOTE: This is to be considered an outdated function and kept in only for reference.
  It's preferable to use cifti_dim_to_nifti_dim() and save_matrix(), which is a _lot_ 
  faster and cleaner.
  """
  
  n_subjects = matrices_array.shape[0]
  n_rows = matrices_array.shape[1]
  n_columns = matrices_array.shape[2]
  
  if coi_indices:
    coi_array = reduce_to_coi(matrices_array, coi_indices)
    if not coi_fname:
      coi_fname = os.path.join(output_path, "coi_matrices.nii.gz")
  
  if clean:
    print("Cleaning matrices before rearranging...")
    clean_matrices_array = np.array(np.empty(matrices_array.shape))
    subject = 0
    while subject < n_subjects:
      matrix_clean = clean_matrix(matrices_array[subject])
      clean_matrices_array[subject] = matrix_clean
      subject += 1
    matrices_array = clean_matrices_array
    if not output_fname:
      output_fname = os.path.join(output_path, "corr_matrices_clean.nii.gz")
  else:
    if not output_fname:
      output_fname = os.path.join(output_path, "corr_matrices_full.nii.gz")
  timer("tic")
  rearray = np.array(np.zeros((n_rows, n_columns, 1, n_subjects)))
  subject = 0
  while subject < n_subjects:
    print("Rearranging matrix of subject {}/{} ({} %).".format(subject+1, n_subjects, round((subject+1)/n_subjects*100, 2)), end="\r")
    row = 0
    while row < n_rows:
      column = 0
      while column < n_columns:
        rearray[row][column][0][subject] = matrices_array[subject][row][column]
        column += 1
      row += 1
    subject += 1 
    sys.stdout.flush # Needed to enable same-line updates on console
  print("\nSaving group NIfTI to {}.".format(output_fname))
  affine = np.array([[-1., 0., 0., n_rows-1], [ 0., 1., 0., -0.], [0., 0., 1., -0.], [0., 0., 0., 1.]])
  save_matrix(rearray, output_fname, affine=affine)
  if coi_fname:
    print("\nCreating COI matrices...")
    coi_array = reduce_to_coi(rearray, coi_indices)
    print("Saving COI group NIfTI to {}.".format(coi_fname))
    save_matrix(coi_array, coi_fname, affine=affine)
  timer("toc", name="rearranging matrices and saving whole-sample NIfTI")

  
def save_matrix(matrix, fname, affine=np.eye(4)):
  image = nib.nifti1.Nifti1Image(matrix, affine=affine)
  nib.save(image, fname)

  
def save_individual_matrices(matrices, subjects, output_dir, clean=False, pconn_dummy=False):
  """
  Take a list of (correlation) matrices, optionally clean them (if we are only interested 
  in the first row and column) and save them as individual NIfTI (and, optionally CIFTI) files
  """
  n_matrices=len(matrices)
  if n_matrices != len(subjects):
    print("ERROR: Mismatch between number of matrices ({}) and number of subjects ({}).".format(n_matrices, len(subjects)))
    exit(1)
  n=0
  if clean:
    clean_matrices = [] # If we want clean matrices, it makes sense to return a list of cleaned matrices
  if pconn_dummy:
    img_dummy = nib.load(pconn_dummy)

  for matrix in matrices:
    fname=os.path.join(output_dir, "{}.nii".format(subjects[n].rstrip('.txt')))
    print("Saving matrix {}/{} to {} ({} %)".format(n+1, n_matrices, fname, round(((n+1)/n_matrices)*100, 2)))
    if clean:
      matrix_to_save = clean_matrix(matrix)
      clean_matrices.append(matrix_clean)
    else:
      matrix_to_save = matrix
    if pconn_dummy:
      img_new_fname = os.path.join(output_dir, "{}.pconn.nii".format(subjects[n]))
      img_new = nib.Cifti2Image(matrix, header=img_dummy.get_header(), file_map=img_dummy.file_map)
      print("\t Saving pconn to {}".format(img_new_fname))
      img_new.to_filename(img_new_fname)
    save_matrix(matrix_to_save, fname)
    n += 1
  if clean:
    return clean_matrices
#  fname = os.path.join(output_path, "{}.nii".format(len(clean_matrices)))
#  print("Saving combined NIfTI file for all matrices (n = {}) to {}.".format(len(clean_matrices), fname))
#  image = nib.cifti2.Cifti2Image(clean_matrices, affine=np.eye(4))
#  nib.save(image, fname)


def compute_correlations(time_series):
  correlation_measure = ConnectivityMeasure(kind='partial correlation')
  print("Fit-transforming time series...")
  timer("tic")
  correlation_matrices = correlation_measure.fit_transform(time_series)
  timer("toc", name="fitting all time series")
  return correlation_measure, correlation_matrices


def clean_matrix(matrix):
  tmp_matrix = matrix.copy()
  # Set all values but first row and first column to zero (our ROI)
  tmp_matrix[1:, 1:] = 0
  # Set correlation of ROI with itself to 0, too
  tmp_matrix[0, 0] = 0 # Actually not needed anymore due to np.fill_diagional later on
  return tmp_matrix


def make_temp_dir():
  tmp_base = "/tmp/roi_connectome"
  if not os.path.exists(tmp_base):
    os.makedirs(tmp_base)
  tmp_dir = mkdtemp(dir=tmp_base)
  return tmp_dir


def plot_all(matrix, title, time_series, coords_file, labels_file, output_dir, clean=True, vmin=None, vmax=None, nan_matrix=False, pconn_dummy=False):
  # Before anything else, save matrix (i.e. corelation_measure.mean) as binary npy 
  # file (e.g. to manually create pconn CIFTIs)
  matrix_fname = os.path.join(output_dir, 'correlation_measure-mean_.npy') 
  np.save(matrix_fname, matrix)
  # Build pconn file
  if pconn_dummy:
    img_new_fname = os.path.join(output_dir, 'correlation_measure-mean_.pconn.nii')
    img_dummy = nib.load(pconn_dummy)
    img_new = nib.Cifti2Image(matrix, header=img_dummy.get_header(), file_map=img_dummy.file_map)
    print("Saving pconn file to {}...".format(img_new_fname))
    img_new.to_filename(img_new_fname)

  plot_title = title + ", n = {}".format(len(time_series))
  plot_title_orig = plot_title # We might modify plot_title later on
  coordinates = np.loadtxt(coords_file)
  labels = load_list_from_file(labels_file)

  if not os.path.exists(output_dir):
    os.makedirs(output_dir)
  
  if nan_matrix is False:
    nan_matrix = matrix
  timer("tic")
  np.fill_diagonal(matrix, 0)
 # plotting.plot_matrix(nan_matrix, colorbar=True, figure=(40, 40), labels=labels, auto_fit=True, vmin=vmin, vmax=vmax)
 # plt.title(plot_title, fontsize=50)
 # fig_fname = "{}/correlation_matrix.svg".format(output_dir)
 # plt.savefig(fig_fname)
 # plt.clf()

  # We need to introduce this as clean_matrix does something unexpected with np.nan values
  # and plot_connectome specifically cannot deal with this
  if clean: 
    print("Cleaning matrix ...")
    mean_roi_matrix = clean_matrix(matrix)
    plot_title = plot_title + " (clean)"
    fname_clean = "clean"
  else:
    mean_roi_matrix = matrix
    fname_clean = "nonclean"
  
  fig_fname = os.path.join(output_dir, "correlation_matrix_{}".format(fname_clean))
  print("Plotting matrix...")
  plotting.plot_matrix(mean_roi_matrix, colorbar=True, figure=(40, 40), labels=labels, auto_fit=True, vmin=vmin, vmax=vmax)
  plt.title(plot_title, fontsize=50)
  print("Saving plot to {}.svg..".format(fig_fname))
  plt.savefig("{}.svg".format(fig_fname))
  #plt.savefig("{}.png".format(fig_fname), dpi=1080) # PNG export does not seem to work in this case. Who knows..
  timer("toc", name="plotting and saving matrices")
  plt.clf()
  print("Plotting connectome ...")
  timer("tic")
  # Manually create new figure because for some reason plotting.plot_connectome() won't
  # accept figure size like plotting.plot_matrix()
  connectome_figure = plt.figure(figsize=(10, 5)) 
  # plot_connectome does not process np.nan values, we still have to pass np.nan_to_num(matrix) to plot_all function
  plotting.plot_connectome(np.nan_to_num(mean_roi_matrix), coordinates, figure=connectome_figure, colorbar=True, 
                           node_size=30, title=plot_title, edge_vmin=vmin, edge_vmax=vmax)
  fig_fname = os.path.join(output_dir, "roi_connectome_{}".format(fname_clean))
  plt.savefig("{}.svg".format(fig_fname))
  plt.savefig("{}.png".format(fig_fname), dpi=1080)
  timer("toc", name="plotting and saving connectome")
  plt.clf()

  html_fname = os.path.join(output_dir, "roi_connectome_90_{}.html".format(fname_clean))
  print("Constructing interactive HTML connectome...")
  timer("tic")
  web_connectome = plotting.view_connectome(mean_roi_matrix, coordinates, edge_threshold="99%", node_size = 6, symmetric_cmap=False)
  html_fname = os.path.join(output_dir, "roi_connectome_{}.html".format(fname_clean))
  web_connectome.save_as_html(html_fname)
  timer("toc", name="plotting and saving HTML connectome")


def clean_palm_results(fname, labels_file, coords_file, alpha=0.95):
  nifti = nib.load(fname)
  data = nifti.get_data()
  labels = load_list_from_file(labels_file)
  coords = load_list_from_file(coords_file)

#  for column_id, column in enumerate(data[0]): # Only the first row is interesting to us
#    if column[0][0] <= alpha: # There is just one "subject" so non-iterating over the last index is fine
#      data[0][column_id][0] = 0 # NA might be better?
#      if column_id > 0:
#        labels[column_id] = ""
  data[data < alpha] = np.nan # Fancy indexing ftw
  nan_indeces = []
  for index, value in enumerate(data[0]): # Get NaN indices to clean labels/coords
    if np.isnan(value):
      if index > 0: # Retain ROI label
        labels[index] = ""
        # Effectively "remove" all NaN nodes by setting their coordinates to our ROI's
        # It works, there has to be a more elegant way, though
        coords[index] = coords[0]
#        coords[index] = "NA"
#  coords = [ coord for coord in coords if coord != "NA" ]
  print("Length of coordinates list: {}".format(len(coords)))
  tmp_dir = make_temp_dir()
  labels_out_fname = os.path.join(tmp_dir, 'clean_labels.txt')
  coords_out_fname = os.path.join(tmp_dir, 'clean_coords.txt')
  with open(labels_out_fname, 'w+') as f:
    for label in labels:
      f.write("%s\n" % label) 
  with open(coords_out_fname, 'w+') as f:
    for label in coords:
      f.write("%s\n" % label)
  # Rearrange NIfTI matrix format to matrix format expected by nilearn
  palm_matrix = np.array(np.zeros((data.shape[3], data.shape[0], data.shape[1])))
  for row_id, row in enumerate(data):
    for column_id, column in enumerate(row):
      for subject_id, subject in enumerate(column):
        palm_matrix[subject_id][row_id][column_id] = subject
        palm_matrix[subject_id][column_id][row_id] = subject # has to be symmetric
  palm_matrix_lite = palm_matrix[~np.isnan(palm_matrix)]
  n_supra_values = len(palm_matrix_lite)
  print("Number of suprathreshold values: {}".format(n_supra_values))
  return palm_matrix, palm_matrix_lite, labels_out_fname, coords_out_fname, n_supra_values


def symmetrize(matrix, mirror_lower=False):
  """
  Return a symmetrized matrix.
  
  If mirror_lower is true, mirror the lower triangle to the upper triangle.
  This method works with matrices containg nan values, too.
  """
  if mirror_lower:
    matrix_symmetric = np.tril(matrix) + np.triu(matrix.T, 1)
  else:
    matrix_symmetric = matrix + matrix.T - np.diag(matrix.diagonal())
  return matrix_symmetric


def nifti_dim_to_cifti_dim(matrix):
  """Converts cols:rows:1:subjects to subjects:cols:rows"""
  matrix_cifti = matrix[:, :, 0, :] # Get rid of dimension #3
  matrix_cifti = np.moveaxis(matrix_cifti, -1, 0) # Put last dimension first
  
  return matrix_cifti


def cifti_dim_to_nifti_dim(matrix):
  """Converts subjects:cols:rows to cols:rows:1:subjects"""
  matrix_nifti = np.moveaxis(matrix, 0, -1) # Subjects:cols:rows -> cols:rows:subjects
  matrix_nifti = np.expand_dims(matrix_nifti, 2) # Cols:rows:subjects -> Colos:rows:1:subjects
  
  return matrix_nifti


def reduce_to_coi(matrices, indices):
    """
    Reduce the number of correlations in a matrices array to correlations of interest (COI).
    This is, for instance, useful to reduce the number of statistical tests we would
    have to correct for later on.
    
    Input: matrices to be "cleaned" and an iterable object of row/column indices to be preserved
    Returns: reduced matrices with all non-interesting values replaced with np.nan
    """
    from matplotlib.cbook import flatten
    
    # Leave the original matrix untouched
    matrices_clean = matrices.copy()
    
    indices_list = [] # Create dummy list to make iteration fool-proof
    indices_list.append(indices)
    indices_flattened = list(flatten(indices_list)) # Unnest nested iterables
    ndim = len(matrices.shape)
    nrow = matrices.shape[1] # Stricly speaking, this would be ncol for NIfTI style
    
    if ndim == 3: # Assume subjects:columns:rows ("CIFTI style")
      print("Assuming CIFTI style with {} rows.".format(nrow))
      for idx in range(0, nrow):
        matrices_clean[:, idx, idx+1:] = np.nan # Nan everything above the diagonal
        if idx not in indices_flattened:
          matrices_clean[:, :, idx] = np.nan # Nan everything except our COIs
    elif ndim == 4: # Assume columns:rows:1:subjects ("NIfTI style")
      print("Assuming NIfTI style with {} cols".format(nrow))
      for idx in range(0, nrow):
        matrices_clean[idx, idx+1:, :, :] = np.nan # Nan everything above the diagonal
        if idx not in indices_flattened:
          matrices_clean[:, idx, :, :] = np.nan # Nan everything else except our COIs      
    else:
      raise ValueError("Number of dimensions of input matrix must be either 3 or 4.")
      
    return matrices_clean
  
def get_nimg_data(fname):
    data = nib.load(fname).get_fdata()
    return data


def wbc_cifti_convert_to_nifti(input_fname, output_fname):
    subprocess.run(["wb_command", "-cifti-convert", "-to-nifti", input_fname, output_fname], capture_output=True)


def run_film_gls(input_fname, output_dir, design_fname, contrast_fname):
    subprocess.run(["film_gls", "--in={}".format(input_fname), "--rn={}".format(os.path.join(output_dir, "film_gls_output")),
                    "--pd={}".format(design_fname), "--con={}".format(contrast_fname), "--thr=1", "--mode=volumetric", "-v"])
    return os.path.join(output_dir, "film_gls_output", "prewhitened_data.nii.gz")

  
def prewhiten_cifti(cifti_fname, design_fname, contrast_fname, tmpdir):
    fakenifti_fname = os.path.join(tmpdir.name, "fakenifti.nii.gz")
    # We could use singlesubject_cifti_to_fake_nifti(), but we have to save a temporary file anyway for film_gls and wb_command is certainly more robust
    wbc_cifti_convert_to_nifti(cifti_fname, fakenifti_fname)
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
                evs = np.loadtxt(ev_fname, ndmin=2) # ndmin=2 added for single-line EV txt files
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
    with open("ts_files", 'w+') as f:
        for filename in fname_list:
            f.write("%s\n" % filename)

            
def singlesubject_cifti_to_fake_nifti(cifti_data):
    #newarray = np.array(np.zeros((cifti_data.shape[1], 1, 1, cifti_data.shape[0])))
    #for idx in range(0, len(cifti_data)): 
    #    for idy in range(0, len(cifti_data[idx])): 
    #        newarray[idy][0][0][idx] = cifti_data[idx][idy]
    newarray = np.transpose(cifti_data) # x:y -> y:x
    newarray = np.expand_dims(newarray, 1) # y:x -> y:1:x
    newarray = np.expand_dims(newarray, 1) # y:1:x -> y:1:1:x
    affine = [[1, 0, 0, 0],[0, 1, 0, 0],[0, 0, 1, 0],[0, 0, 0, 1]]
    img = nib.nifti1.Nifti1Image(newarray, affine=affine)
    
    return img

  
def plot_palm_new(palm_results_fname, title, coords, labels, alpha=1.3, scale=False):
    data = get_nimg_data(palm_results_fname)
    adjmatrix = data[:, :, 0, 0]
    # Plot all p values
    fig_matrix = plotting.plot_matrix(adjmatrix, colorbar=True, figure=(40, 40), labels=labels, auto_fit=True)
    fig_connectome = plotting.plot_connectome(symmetrize(adjmatrix, mirror_lower=True), coords, colorbar=True, node_size=15, title=title)
    web_connectome = plotting.view_connectome(adjmatrix, coords, node_size = 6, symmetric_cmap=False)
    
    # Check if there are any suprathreshold values
    nsupthr = len(adjmatrix[adjmatrix >= alpha])
    if nsupthr == 0:
      print("\nNo values survive the threshold of {}.".format(alpha))
      fig_matrix_clean, fig_connectome_clean, web_connectome_clean = None, None, None
    else:
      print("Number of values to survive the threshold of {}: {}".format(alpha, nsupthr))  
      # Purge matrix, coordinates and labels of subthreshold values
      adjmatrix[adjmatrix < alpha] = np.nan
      supthr_indices = np.argwhere(~np.isnan(adjmatrix))
      labels_clean = ["" if i not in supthr_indices else x for i, x in enumerate(labels)]
      coords_clean = [[np.nan, np.nan, np.nan] if i not in supthr_indices else x for i, x in enumerate(coords)]

      fig_matrix_clean = plotting.plot_matrix(adjmatrix, colorbar=True, figure=(40, 40), labels=labels_clean, auto_fit=True)
      fig_connectome_clean = plotting.plot_connectome(symmetrize(np.nan_to_num(adjmatrix), mirror_lower=True), coords_clean, figure=plt.figure(figsize=(10, 5)), colorbar=True, node_size=30, title=title)
      web_connectome_clean = plotting.view_connectome(adjmatrix, coords_clean, node_size = 6, symmetric_cmap=False)
    
    return fig_matrix, fig_matrix_clean, fig_connectome, fig_connectome_clean, web_connectome, web_connectome_clean