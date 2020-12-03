#!/usr/bin/python3
# Adapted from https://github.com/esfinn/cpm_tutorial/blob/master/cpm_tutorial.ipynb

import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
from pathlib import Path
import os
import warnings
import nibabel as nib
import sys
from nilearn import plotting
sys.path.append('/home/test/Projekte/diss/hcp-suite/lib')
from hcpsuite import get_nimg_data, printv, nifti_dim_to_cifti_dim

global verbose # We need this as global variables are actually module-level, not truly global
verbose = True

def convert_matrices_to_dataframe(array, subj_ids):
  """
  Takes a NumPy array (subjects:parcels:parcels)and converts into a Pandas dataframe fit 
  for downstream CPM analyses
  """
  assert array.shape[0] == len(subj_ids), "Number of subject IDs is not equal to number of subjects in neuroimage file"
  all_fc_data = {}
  n = 0
  for id in subj_ids:
    printv("Flattening matrix of subject {} ({} of {}...)".format(id, n+1, len(subj_ids)), update=True)
    tmp = array[n, :, :] # Get matrix of a single subject
    all_fc_data[id] = tmp[np.triu_indices_from(tmp, k=1)] # Only use upper triangle of symmetric matrix
    n += 1
  printv("\nCreating DataFrame from matrices...")
  all_fc_data = pd.DataFrame.from_dict(all_fc_data, orient='index')
  
  return all_fc_data


def get_matrices_from_nimg_file(nimg_file, subj_ids):
  """
  Reads a neuroimage file (e.g. in NIfTI oder CIFTI format) containing connectivity
  matrices in the conventional shape of parcels:parcels:1:subjects and returns a dataframe
  of flattened matrices, i.e. functional connectivity per subject (each row represents
  functional connectivity of one subject, first column is subject ID).
  
  The neuroimage file's order of subjects _must_ be the same as the subject IDs list's
  """
 
  matrices = get_nimg_data(nimg_file)[:, :, 0, :] # -> parcels:parcels:subjects
  assert matrices.shape[2] == len(subj_ids), "Number of subject IDs is not equal to number of subjects in neuroimage file"
  matrices_reshaped = np.transpose(matrices, (2, 0, 1)) # -> subjects:parcels:parcels
  all_fc_data = convert_matrices_to_dataframe(matrices_reshaped, subj_ids)
  
  return all_fc_data


def get_behav_data(fname, ids):
  """
  Read behavioural data for specified IDs from CSV and return dataframe
  """
  
  all_behav_data = pd.read_csv(fname, dtype={'Subject': str}, low_memory = False)
  all_behav_data.set_index('Subject', inplace=True)
  all_behav_data = all_behav_data.loc[ids] # Keep only data from specified subjects
  
  return all_behav_data


def create_clean_upper_reshaped(matrix):
  """
  Takes a matrix in the shape of parcels:parcels:1:subjects (e.g. NIfTI style)
  and converts it into a matrix of subjects:parcels:parcels, and cleans it of
  edges of non-interest
  """
  matrix = nifti_dim_to_cifti_dim(matrix)
  matrix[:, 99:512, 99:512] = 0
  matrix = np.triu(matrix)
  
  return matrix
  

def create_kfold_indices(subj_list, k = 10):
  """
  Splits list of subjects into k folds for cross-validation.
  """
  
  n_subs = len(subj_list)
  n_subs_per_fold = n_subs//k # floor integer for n_subs_per_fold

  indices = [[fold_no]*n_subs_per_fold for fold_no in range(k)] # generate repmat list of indices
  remainder = n_subs % k # figure out how many subs are left over
  remainder_inds = list(range(remainder))
  indices = [item for sublist in indices for item in sublist]    
  [indices.append(ind) for ind in remainder_inds] # add indices for remainder subs

  assert len(indices)==n_subs, "Length of indices list does not equal number of subjects"

  np.random.shuffle(indices) # shuffles in place

  return np.array(indices)

  
def split_train_test(subj_list, indices, test_fold):
  """
  For a subj list, k-fold indices, and given fold, returns lists of train_subs and test_subs
  """

  train_inds = np.where(indices!=test_fold)
  test_inds = np.where(indices==test_fold)

  train_subs = []
  for sub in subj_list[train_inds]:
    train_subs.append(sub)

  test_subs = []
  for sub in subj_list[test_inds]:
    test_subs.append(sub)

  return (train_subs, test_subs)


def get_train_test_data(all_fc_data, train_subs, test_subs, behav_data, behav):
  """
  Extracts requested FC and behavioral data for a list of train_subs and test_subs
  """

  train_vcts = all_fc_data.loc[train_subs, :]
  test_vcts = all_fc_data.loc[test_subs, :]

  train_behav = behav_data.loc[train_subs, behav]

  return (train_vcts, train_behav, test_vcts)


def select_features(train_vcts, train_behav, r_thresh=0.2, corr_type='pearson'):
  """
  Runs the CPM feature selection step: 
  - correlates each edge with behavior, and returns a mask of edges that are correlated above some threshold, one for each tail (positive and negative)
  """
  global verbose
  assert train_vcts.index.equals(train_behav.index), "Row indices of FC vcts and behavior don't match!"

  # Correlate all edges with behav vector
  if corr_type =='pearson':
    cov = np.dot(train_behav.T - train_behav.mean(), train_vcts - train_vcts.mean(axis=0)) / (train_behav.shape[0]-1)
    corr = cov / np.sqrt(np.var(train_behav, ddof=1) * np.var(train_vcts, axis=0, ddof=1))
  elif corr_type =='spearman':
    corr = []
    for edge in train_vcts.columns:
      r_val = sp.stats.spearmanr(train_vcts.loc[:,edge], train_behav)[0]
      corr.append(r_val)

  # Define positive and negative masks
  mask_dict = {}
  mask_dict["pos"] = corr > r_thresh
  mask_dict["neg"] = corr < -r_thresh
    
  printv("  - Found ({}/{}) edges positively/negatively correlated (threshold: {}) with behavior in the training set".format(mask_dict["pos"].sum(), mask_dict["neg"].sum(), r_thresh)) # for debugging
  printv("  - Max r pos: {}, max r neg: {}".format(corr.max(), corr.min()))

  return mask_dict


def is_invalid_array(array, sanitize = False):
  """Takes an array and checks the following:
       - all zeros?
       - any infinity values?
       - any NaN values?
     If any of these conditions are met, it returns true and a message indicating the fulfilled condition.
     If sanitize = True, also return sanitized array
  """
  
  is_invalid = False
  msg = "Array valid."
  
  if np.isnan(array).any():
    is_invalid = True
    msg = "The array contains NaN values"
    #if sanitize: Do we want this?
    #  array[np.isnan(array)] = 0
  elif np.isinf(array).any():
    is_invalid = True
    msg = "The array contains infinite values"
    #if sanitize: # Do we want this?
    #  array[np.isposinf(array)] = 1
    #  array[np.isneginf(array)] = -1
  elif not np.any(array):
    is_invalid = True
    msg = "The array contains all zeros"
    if sanitize:
      array[array.all()] = np.nan
  
  if sanitize:
    return is_invalid, msg, array
  else:
    return is_invalid, msg
    

def build_model(train_vcts, mask_dict, train_behav):
  """
  Takes a feature mask, sums all edges in the mask for each subject, and uses simple linear 
  regression to relate summed network strength to behavior; returns a dictionary with the model
  """

  assert train_vcts.index.equals(train_behav.index), "Row indices of FC vcts and behavior don't match!"

  model_dict = {}

  # Loop through pos and neg tails
  X_glm = np.zeros((train_vcts.shape[0], len(mask_dict.items())))

  t = 0
  for tail, mask in mask_dict.items():
    X = train_vcts.values[:, mask].sum(axis=1)
    X_glm[:, t] = X
    y = train_behav
    is_invalid, invalid_msg = is_invalid_array(X)
    if is_invalid:
      print("X: {}".format(invalid_msg))
      return False
    is_invalid, invalid_msg = is_invalid_array(y)
    if is_invalid:
      print("y: {}".format(invalid_msg))
      return False
    
    (slope, intercept) = np.polyfit(X, y, 1)
    model_dict[tail] = (slope, intercept)
    t += 1
        
  X_glm = np.c_[X_glm, np.ones(X_glm.shape[0])]
  model_dict["glm"] = tuple(np.linalg.lstsq(X_glm, y, rcond=None)[0])

  return model_dict

def apply_model(test_vcts, mask_dict, model_dict):
  """
  Applies a previously trained linear regression model to a test set to generate predictions of behavior.
  """

  behav_pred = {}

  X_glm = np.zeros((test_vcts.shape[0], len(mask_dict.items())))

  # Loop through pos and neg tails
  t = 0
  for tail, mask in mask_dict.items():
    X = test_vcts.loc[:, mask].sum(axis=1)
    X_glm[:, t] = X

    slope, intercept = model_dict[tail]
    behav_pred[tail] = slope*X + intercept
    t+=1

  X_glm = np.c_[X_glm, np.ones(X_glm.shape[0])]
  behav_pred["glm"] = np.dot(X_glm, model_dict["glm"])

  return behav_pred

def plot_predictions(behav_obs_pred, tail="glm"):
    x = behav_obs_pred.filter(regex=("obs")).astype(float)
    y = behav_obs_pred.filter(regex=(tail)).astype(float)

    g = sns.regplot(x=x.T.squeeze(), y=y.T.squeeze(), color='gray')
    ax_min = min(min(g.get_xlim()), min(g.get_ylim()))
    ax_max = max(max(g.get_xlim()), max(g.get_ylim()))
    g.set_xlim(ax_min, ax_max)
    g.set_ylim(ax_min, ax_max)
    g.set_aspect('equal', adjustable='box')
    
    r = get_r_value(behav_obs_pred)
    g.annotate('r = {0:.2f}'.format(r[0]), xy = (0.7, 0.1), xycoords = 'axes fraction')
    
    return g
  
def get_r_value(behav_obs_pred, tail="glm"):
  x = behav_obs_pred.filter(regex=("obs")).astype(float)
  y = behav_obs_pred.filter(regex=(tail)).astype(float)
  
  x[np.isnan(x)] = 0
  y[np.isnan(y)] = 0
  
  r = sp.stats.pearsonr(x.iloc[:, 0], y.iloc[:, 0])
  
  return r
  
  
def plot_consistent_edges(all_masks, tail, thresh = 1., color='gray', coords=None):
    edge_frac = (all_masks[tail].sum(axis=0))/(all_masks[tail].shape[0])
    print("For the {} tail, {} edges were selected in at least {}% of folds".format(tail, (edge_frac>=thresh).sum(), thresh*100))
    edge_frac_square = sp.spatial.distance.squareform(edge_frac)

    node_mask = np.amax(edge_frac_square, axis=0) >= thresh # find nodes that have at least one edge that passes the threshold
    node_size = edge_frac_square.sum(axis=0)*node_mask*20 # size nodes based on how many suprathreshold edges they have

    plotting.plot_connectome(adjacency_matrix=edge_frac_square, edge_threshold=thresh,
                    node_color = color,
                    node_coords=coords, node_size=node_size,
                    display_mode= 'lzry',
                    edge_kwargs={"linewidth": 1, 'color': color})

def plot_consistent_edges_loo(posmasks, thresh=0.13, consistency=0.8, coords=None):
  """Plot edges obtained in a leave-one out CPM above a defined threshold that 
  are selected in at least a defined percentage of subjects."""
  
  posmasks = symmetrize_matrices(posmasks, mirror_lower=True) # Symmetrize matrices
  posmasks_flat = posmasks.reshape(posmasks.shape[0], posmasks.shape[1]*posmasks.shape[2]) # Flatten matrix
  posmasks_flat
  edges_count = np.zeros(posmasks_flat.shape[1])
  edges_mask = np.zeros(posmasks_flat.shape[1])
  for i in range(0, posmasks_flat.shape[0]):
    edges_count[posmasks_flat[i] > thresh] += 1 # Count number of times an edge is suprathreshold
  
  edges_mask[edges_count > posmasks.shape[0]*consistency] = 1
  nodes_mask = edges_mask.reshape((posmasks.shape[1], posmasks.shape[2]))
  printv("There are {} suprathreshold (> {}) edges in {} % of the subjects".format(edges_mask[edges_mask ==1].sum()/2, thresh, consistency*100))
  
  degree_list = []
  for node in range(513):
    degree = nodes_mask[node, :].sum() # Determine degre of each node and add it to list
    degree_list.append(degree)
 
  plotting.plot_connectome(nodes_mask, node_coords=coords, display_mode='lzry', node_size=[degree*20 for degree in degree_list], edge_kwargs={"linewidth": 2})

  
def perform_cpm(all_fc_data, all_behav_data, behav, k=10, **cpm_kwargs):
  """
  Takes functional connectivity and behaviour dataframes, selects a behaviour
  """
  assert all_fc_data.index.equals(all_behav_data.index), "Row (subject) indices of FC vcts and behavior don't match!"

  subj_list = all_fc_data.index # get subj_list from df index
    
  indices = create_kfold_indices(subj_list, k=k)
    
  # Initialize df for storing observed and predicted behavior
  col_list = []
  for tail in ["pos", "neg", "glm"]:
    col_list.append(behav + " predicted (" + tail + ")")
  col_list.append(behav + " observed")
  behav_obs_pred = pd.DataFrame(index=subj_list, columns = col_list)
    
  # Initialize array for storing feature masks
  n_edges = all_fc_data.shape[1]
  all_masks = {}
  all_masks["pos"] = np.zeros((k, n_edges))
  all_masks["neg"] = np.zeros((k, n_edges))
  
  n_folds_completed = 0
  for fold in range(k):
    print("Doing fold {} of {} (successful folds: {})...".format(fold + 1, k, n_folds_completed))
    train_subs, test_subs = split_train_test(subj_list, indices, test_fold=fold)
    train_vcts, train_behav, test_vcts = get_train_test_data(all_fc_data, train_subs, test_subs, all_behav_data, behav=behav)
    mask_dict = select_features(train_vcts, train_behav, **cpm_kwargs)
    all_masks["pos"][fold,:] = mask_dict["pos"]
    all_masks["neg"][fold,:] = mask_dict["neg"]
    model_dict = build_model(train_vcts, mask_dict, train_behav)
    if not model_dict: # build_model returns False instead of a dict if an array is not valid
      print("  - Fold failed -> continuing with next fold...")
      continue # Skip fold if generated arrays are not valid
    behav_pred = apply_model(test_vcts, mask_dict, model_dict)
    for tail, predictions in behav_pred.items():
      behav_obs_pred.loc[test_subs, behav + " predicted (" + tail + ")"] = predictions
    n_folds_completed += 1
  
  print("\nCPM completed. Successful folds: {}".format(n_folds_completed))
  behav_obs_pred.loc[subj_list, behav + " observed"] = all_behav_data[behav]
    
  return behav_obs_pred, all_masks

