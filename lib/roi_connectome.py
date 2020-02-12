#!/usr/bin/python3

import sys
import os
from tempfile import mkdtemp
import time
import numpy as np
import nibabel as nib
from nilearn.connectome import ConnectivityMeasure
from nilearn import plotting
from matplotlib import pyplot as plt
from pandas import read_csv

def usage(argv):
  if len(argv) < 5:
    print("\n\n  Usage: {} <list of time series files> <coordinates file> <labels file> <title> [<output dir>]\n".format(os.path.basename(argv[0])))
    exit(1)
  else:
    ts_file = argv[1]
    coords_file = argv[2]
    labels_file = argv[3]
    title = os.path.basename(argv[4])
    if len(argv) == 6:
      output_dir = argv[5]
      if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    elif len(argv) == 5:
      output_dir = make_temp_dir()
    return ts_file, coords_file, labels_file, title, output_dir

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

def save_whole_sample_nifti(matrices_array, output_path, clean=True):
  """Take an array of matrices (e.g. output from ConnecitivityMeasure()),
  rearrange it to build a multi-subject NIfTI, e.g. for 361 ROIs and 999 subjects:
  (999, 361, 361) -> (361, 361, 1, 999)"""
  n_subjects = matrices_array.shape[0]
  n_rows = matrices_array.shape[1]
  n_columns = matrices_array.shape[2]
  if clean:
    print("Cleaning matrices before rearranging...")
    clean_matrices_array = np.array(np.empty(matrices_array.shape))
    subject = 0
    while subject < n_subjects:
      matrix_clean = clean_matrix(matrices_array[subject])
      clean_matrices_array[subject] = matrix_clean
      subject += 1
    matrices_array = clean_matrices_array
    fname = os.path.join(output_path, "corr_matrices_clean.nii.gz")
  else:
    fname = os.path.join(output_path, "corr_matrices_full.nii.gz")
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
  print("\nSaving group NIfTI to {}.".format(fname))
  affine = np.array([[-1., 0., 0., n_rows-1], [ 0., 1., 0., -0.], [0., 0., 1., -0.], [0., 0., 0., 1.]])
  save_matrix(rearray, fname, affine=affine)
  timer("toc", name="rearranging matrices and saving whole-sample NIfTI")

def save_matrix(matrix, fname, affine=np.eye(4)):
  image = nib.nifti1.Nifti1Image(matrix, affine=affine)
  nib.save(image, fname)

def save_individual_matrices(matrices, subjects, output_path):
  """Take a list of (correlation) matrices, clean them (we are only interested in the first 
  row and column) and save them as individual NIfTI files"""
  n_matrices=len(matrices)
  if n_matrices != len(subjects):
    print("ERROR: Mismatch between number of matrices ({}) and number of subjects ({}).".format(n_matrices, len(subjects)))
    exit(1)
  n=0
  clean_matrices = []
  for matrix in matrices:
    fname=os.path.join(output_path, "{}.nii".format(subjects[n].rstrip('.txt')))
    print("Saving matrix {}/{} to {} ({} %)".format(n+1, n_matrices, fname, round(((n+1)/n_matrices)*100, 2)))
    matrix_clean = clean_matrix(matrix)
    clean_matrices.append(matrix_clean)
    save_matrix(matrix_clean, fname)
    n += 1
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
    img_new.to_filename(img_new_fname)

  plot_title = title + ", n = {}".format(len(time_series))
  plot_title_orig = plot_title # We might modify plot_title later on
  coordinates = np.loadtxt(coords_file)
  labels = load_list_from_file(labels_file)

  if not os.path.exists(output_dir):
    os.makedirs(output_dir)
  
  if nan_matrix is False:
    nan_matrix = matrix
  print("Plotting matrices ...")
  timer("tic")
  np.fill_diagonal(matrix, 0)
  plotting.plot_matrix(nan_matrix, colorbar=True, figure=(40, 40), labels=labels, auto_fit=True, vmin=vmin, vmax=vmax)
  plt.title(plot_title, fontsize=50)
  fig_fname = "{}/correlation_matrix.svg".format(output_dir)
  plt.savefig(fig_fname)
  plt.clf()

  # We need to introduce this as clean_matrix does something unexpected with np.nan values
  # and plot_connectome specifically cannot deal with this
  if clean: 
    print("Cleaning matrix ...")
    mean_roi_matrix = clean_matrix(matrix)
    plot_title = plot_title + " (clean)"
  else:
    mean_roi_matrix = matrix

  plotting.plot_matrix(mean_roi_matrix, colorbar=True, figure=(40, 40), labels=labels, auto_fit=True, vmin=vmin, vmax=vmax)
  plt.title(plot_title, fontsize=50)
  fig_fname = "{}/correlation_matrix_clean.svg".format(output_dir)
  plt.savefig(fig_fname)
  timer("toc", name="plotting and saving matrices")
  plt.clf()
  print("Plotting connectome ...")
  timer("tic")
  # Manually create new figure because for some reason plotting.plot_connectome() won't
  # accept figure size like plotting.plot_matrix()
  connectome_figure = plt.figure(figsize=(10, 5)) 
  # plot_connectome does not process np.nan values, we still have to pass np.nan_to_num(matrix) to plot_all function
  plotting.plot_connectome(np.nan_to_num(mean_roi_matrix), coordinates, figure=connectome_figure, colorbar=True, node_size=30, title=plot_title, edge_vmin=vmin, edge_vmax=vmax)
  fig_fname = "{}/roi_connectome.svg".format(output_dir)
  plt.savefig(fig_fname)
  timer("toc", name="plotting and saving connectome")
  plt.clf()
  
  # Also plot non-clean matrix
  connectome_figure = plt.figure(figsize=(10, 5)) 
  plotting.plot_connectome(np.nan_to_num(matrix), coordinates, figure=connectome_figure, colorbar=True, node_size=30, title=plot_title_orig, edge_vmin=vmin, edge_vmax=vmax)
  fig_fname = "{}/roi_connectome_nonclean.svg".format(output_dir)
  plt.savefig(fig_fname)
  # Non-clean SVGs tend to be ressource-hungry; PNG might be more practical
  fig_fname = "{}/roi_connectome_nonclean.png".format(output_dir)
  plt.savefig(fig_fname)
  timer("toc", name="plotting and saving connectome")
  plt.clf()

  html_fname = "{}/roi_connectome_90.html".format(output_dir)
  print("Constructing interactive HTML connectome...")
  timer("tic")
  web_connectome = plotting.view_connectome(mean_roi_matrix, coordinates, edge_threshold="99%", node_size = 6, symmetric_cmap=False)
  html_fname = "{}/roi_connectome.html".format(output_dir)
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

def symmetrize(matrix):
  return matrix + matrix.T - np.diag(matrix.diagonal())
