import sys
import ray
from ray.util.queue import Queue
import numpy as np
import pandas as pd
from time import sleep
import socket
import os
import pickle
import seaborn as sns
import scipy as sp
from sklearn.model_selection import KFold
from pingouin import partial_corr
from sklearn.linear_model import LinearRegression
from nilearn import plotting
from hcpsuite import load_list_from_file, get_nimg_data, printv, nifti_dim_to_cifti_dim, symmetrize_matrices, symmetrize, timer

def get_behav_data(fname, ids=None):
  """
  Read behavioural data for specified IDs from CSV and return dataframe
  """
  
  behav_data = pd.read_csv(fname, dtype={'Subject': str}, low_memory = False)
  behav_data.set_index('Subject', inplace=True)
  if ids:
    behav_data = behav_data.loc[ids] # Keep only data from specified subjects
  
  return behav_data

def convert_matrices_to_dataframe(array, subj_ids):
  """
  Takes a NumPy array (subjects:parcels:parcels) and converts it into a Pandas dataframe fit 
  for downstream CPM analyses
  """
  assert array.shape[0] == len(subj_ids), "Number of subject IDs is not equal to number of subjects in neuroimage file"
  fc_data = {}
  n = 0
  for id in subj_ids:
    printv("Flattening matrix of subject {} ({} of {}...)".format(id, n+1, len(subj_ids)), update=True)
    tmp = array[n, :, :] # Get matrix of a single subject
    fc_data[id] = tmp[np.triu_indices_from(tmp, k=1)] # Only use upper triangle of symmetric matrix
    n += 1
  printv("\nCreating DataFrame from matrices...")
  fc_data = pd.DataFrame.from_dict(fc_data, orient='index')
  
  return fc_data

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
  

def plot_consistent_edges(all_masks, tail, thresh = 1., color='gray', coords=None, **plot_kwargs):
    """Plots edges which are consistent in a defined percentage of folds"""
    edge_frac = (all_masks[tail].sum(axis=0))/(all_masks[tail].shape[0])
    summary = "{} suprathreshold (>= {} %) edges in {} % of the subjects".format((edge_frac >= thresh).sum(), thresh*100, "TODO")
    print("For the {} tail, {} edges were selected in at least {}% of folds".format(tail, (edge_frac>=thresh).sum(), thresh*100))
    edge_frac_square = sp.spatial.distance.squareform(edge_frac)

    node_mask = np.amax(edge_frac_square, axis=0) >= thresh # find nodes that have at least one edge that passes the threshold
    node_size = edge_frac_square.sum(axis=0)*node_mask*20 # size nodes based on how many suprathreshold edges they have

    plotting.plot_connectome(adjacency_matrix=edge_frac_square, edge_threshold=thresh,
                    node_color = color,
                    node_coords=coords, node_size=node_size,
                    display_mode= 'lzry',
                    edge_kwargs={"linewidth": 1, 'color': color},
                    **plot_kwargs)
    
    return (edge_frac_square >= thresh).astype(int), node_mask

def plot_top_n_edges(matrix, top_n, img_base=None, img_naming_scheme='label', img_suffix='png', labels=None, color='red', check_for_symmetry=False):
    """
    Takes a connectivity/correlation/prediction matrix (e.g. a meaned and masked r_mat from
    pcorr_dfs_to_rmat_pmat(), and plots the top n edges as NetworkX graphs. img_naming_scheme
    is either 'label' (in which case a list of labels has to be supplied) or 'number'
    """
    import networkx as nx
    from matplotlib import pyplot as plt
    
    if img_naming_scheme == 'label' and not labels:
        raise ValueError("A list of labels needs to be specified if img_namingscheme is 'label'")
    elif labels:
        assert len(labels) == matrix.shape[0], "Number of labels does not match number of nodes"
    
    if check_for_symmetry: # is_symmetric() is somewhat unreliable, therefore it is disabled by default
        if is_symmetric(matrix):
            matrix = np.tril(matrix) # Clear upper triangle in symmetric matrix to avoid duplicate edges
    else:
            matrix = np.tril(matrix) # Clear upper triangle in symmetric matrix to avoid duplicate edges
    matrix = abs(matrix) # Let this function work for both pos and neg tails
    top_n_indices = np.unravel_index(np.argpartition(-matrix.flatten(), range(top_n))[:top_n], matrix.shape)
    top_n_values = matrix[top_n_indices]
    figures = []
    for n in range(top_n):
        G = nx.Graph()
        nodes = []
        nodes.append(top_n_indices[0][n])
        nodes.append(top_n_indices[1][n])
        for node in nodes:
            G.add_node(node)
            if labels:
                G.nodes[node]['label'] = labels[node]
            if img_base:
                if img_naming_scheme == 'label':
                    fname = os.path.join(img_base, "{}.{}".format(labels[node], img_suffix))
                elif img_naming_scheme == 'number':
                    fname = os.path.join(img_base, "{}.{}".format(node+1, img_suffix)) # NumPy arrays are 0-indexed, nodes are not
                else:
                    raise ValueError("img_naming_scheme must be either 'label' or 'number'")
                img = plt.imread(fname)
                G.nodes[node]['image'] = img
        G.add_edge(nodes[0], nodes[1], weight=round(top_n_values[n],6))
        
        # The next section relies heavily on https://stackoverflow.com/a/53968787      
        pos = nx.planar_layout(G)

        fig = plt.figure(figsize=(8,4))
        ax = plt.subplot(111)
        ax.set_aspect('equal')
        edge_labels = nx.get_edge_attributes(G, 'weight')
        nx.draw_networkx_edges(G, pos, ax=ax, edge_color=color)
        nx.draw_networkx_edge_labels(G, pos, ax=ax, edge_labels=edge_labels) 

        plt.xlim(-0.9,0.9)
        plt.ylim(0.5,-0.5)

        trans = ax.transData.transform
        trans2 = fig.transFigure.inverted().transform

        img_size = 0.5 # this is the image size
        img_size_12 =img_size/2.0
        for n in G:
            xx, yy = trans(pos[n]) # figure coordinates
            xa, ya = trans2((xx,yy)) # axes coordinates
            a = plt.axes([xa-img_size_12,ya-img_size_12, img_size, img_size])
            a.set_aspect('equal')
            a.imshow(G.nodes[n]['image'])
            a.text(0.5, -0.2, G.nodes[n]['label'], transform=a.transAxes, horizontalalignment='center')
            a.axis('off')
        #ax.text(0, 0.1, G.)
        ax.axis('off')
        figures.append(fig)
        
    return figures

def plot_top_n_nodes(degrees, top_n, parcellation_img, bg_img, n_slices=5, display_modes=['x', 'y', 'z']):
    n_nodes = len(degrees)
    top_n_nodes = np.argsort(degrees)[-top_n:]
    
    new_data = np.zeros(parcellation_img.shape)
    for node in top_n_nodes:
        # In a parcellation image, parcels = nodes are represented as neighbouring 
        # voxels with values corresponding to the parcel/node number
        new_data[rtg_data == node] = node
    new_img = nib.Nifti1Image(new_data, rtg.affine, rtg.header)
    
    figures_dict = {}
    for display_mode in display_modes:
        plotting.plot_roi(new_img, cut_coords=5, display_mode=display_mode, bg_img=bg_img,     black_bg=False)
        figures_dict[display_mode] = plt.gcf()
        
    return figures_dict
    

def plot_predictions(predictions, tail="glm", save=False, title=None, fname='predictions.svg', color='gray'):
    x = predictions['observed'].astype(float)
    y = predictions[tail].astype(float)

    g = sns.regplot(x=x.T.squeeze(), y=y.T.squeeze(), color=color)
    ax_min = min(min(g.get_xlim()), min(g.get_ylim()))
    ax_max = max(max(g.get_xlim()), max(g.get_ylim()))
    g.set_xlim(ax_min, ax_max)
    g.set_ylim(ax_min, ax_max)
    g.set_aspect('equal', adjustable='box')
    g.set_title(title)
    
    r = get_r_value(x, y)
    g.annotate('r = {0:.2f}'.format(r), xy = (0.7, 0.1), xycoords = 'axes fraction')
    
    if save:
      fig = g.get_figure()
      fig.savefig(fname, bbox_inches="tight")
    
    return g
  
def get_r_value(x, y):
  x[np.isnan(x)] = 0
  y[np.isnan(y)] = 0
  
  r = sp.stats.pearsonr(x, y)[0]
  
  return r

def get_kfold_indices(subs_list, k):
  """
  Returns a dictionary of train and test indices of length k
  """
  subs_array = np.array(subs_list)
  kf = KFold(n_splits=k, shuffle=True)
  kfold_indices = {}
  kfold_indices['train'] = []
  kfold_indices['test'] = []
  for train_index, test_index in kf.split(subs_list):
    kfold_indices['train'].append(subs_array[train_index]) 
    kfold_indices['test'].append(subs_array[test_index])
  
  return kfold_indices

def clean_kfold_indices(kfold_indices, behav_data, noneb_group='Mother_ID'):
  """
  Cleans the train part of kfold indices of entries sharing the same group.
  This was written to exclude closely related subjects in HCP data by using 
  'Mother_ID' as the noneb_group (read: non-exchangeability group).
  
  Returns the cleaned kfold_indices dictionary.
  """
  kfold_indices_clean = {} # Rebuilding from scratch seems easier than using copy.deepcopy()
  kfold_indices_clean['train'] = []
  kfold_indices_clean['test'] = []
  for fold in range(len(kfold_indices['train'])):
    noneb_values = behav_data.loc[kfold_indices['test'][fold], noneb_group]
    mask = behav_data.loc[kfold_indices['train'][fold], noneb_group].isin(noneb_values)
    kfold_indices_clean['train'].append(kfold_indices['train'][fold][~mask])
    kfold_indices_clean['test'].append(kfold_indices['test'][fold])

  return kfold_indices_clean


def pcorr_dfs_to_rmat_pmat(pcorr_dfs):
  """
  Takes a list of correlation DataFrames (post-feature selection), extracts their
  p and r values and puts them in corresponding matrices.
  
  Input: List of DataFrames
  Returns: r_mat and p_mat
  """
  n_mats = len(pcorr_dfs) # Number of individual matrices
  n_edges = len(pcorr_dfs[0]['r'])
  shape_mats = sp.spatial.distance.squareform(np.zeros(n_edges)).shape # Get squared dimension of matrix
  r_mat = np.zeros((n_mats, shape_mats[0], shape_mats[1]))
  p_mat = np.zeros((n_mats, shape_mats[0], shape_mats[1]))

  for n in range(n_mats):
    r_mat[n] = sp.spatial.distance.squareform(pcorr_dfs[n]['r'])
    p_mat[n] = sp.spatial.distance.squareform(pcorr_dfs[n]['p-val'])
  
  return r_mat, p_mat
  

def get_suprathr_edges_new(df_dict, p_thresh_pos=None, p_thresh_neg=None, r_thresh_pos=None, r_thresh_neg=None, percentile_neg=None, percentile_pos=None, top_n_pos=None, top_n_neg=None):
  folds_list = list(df_dict.keys())
  n_edges = len(df_dict[folds_list[0]])
  masks_dict = {}

  for fold in folds_list:
    pcorr_df = df_dict[fold]
    n_edges = len(df_dict[fold])
    masks_dict[fold] = {}
    suprathr_edges_mask = {}
    if p_thresh_pos and p_thresh_neg:
      suprathr_edges_mask['pos'] = (pcorr_df['r'] > 0) & (pcorr_df['p-val'] <= p_thresh_pos)
      suprathr_edges_mask['neg'] = (pcorr_df['r'] < 0) & (pcorr_df['p-val'] <= p_thresh_neg)
    elif r_thresh_pos and r_thresh_neg:
      suprathr_edges_mask['pos'] = pcorr_df['r'] > r_thresh_pos
      suprathr_edges_mask['neg'] = pcorr_df['r'] < -abs(r_thresh_neg) # r_thresh_neg can be both given as a positive or a negative value
    elif percentile_pos and percentile_neg:
      r_thresh_pos = np.nanpercentile(pcorr_df['r'], percentile_pos)
      r_thresh_neg = np.nanpercentile(pcorr_df['r'][pcorr_df['r'] < 0], 100 - percentile_neg)
      suprathr_edges_mask['pos'] = pcorr_df['r'] > r_thresh_pos
      suprathr_edges_mask['neg'] = pcorr_df['r'] < -abs(r_thresh_neg)
    elif top_n_pos and top_n_neg:
      suprathr_edges_mask['pos'] = np.zeros(pcorr_df.shape[0])
      suprathr_edges_mask['neg'] = np.zeros(pcorr_df.shape[0])
      suprathr_edges_mask['pos'][np.argpartition(pcorr_df['r'][pcorr_df['r'].notna()], -top_n_pos)[-top_n_pos:]] = 1
      suprathr_edges_mask['neg'][np.argpartition(pcorr_df['r'][pcorr_df['r'].notna()], top_n_neg)[:top_n_neg]] = 1
    else:
      raise TypeError('Either p_thresh_{neg, pos} or r_thresh_{neg, pos} or percentile_{neg, pos} or top_n_{pos, neg} needs to be defined.')

    printv("Fold {}: Pos/neg suprathreshold edges (max r pos/max r neg): {}/{} ({}/{})".format(fold+1, suprathr_edges_mask['pos'].sum(), suprathr_edges_mask['neg'].sum(), pcorr_df['r'].max(), pcorr_df['r'].min()))
    for tail in ('pos', 'neg'):
        masks_dict[fold][tail] = np.zeros(n_edges)
        masks_dict[fold][tail][:] = suprathr_edges_mask[tail].astype(bool)

  return masks_dict

def get_suprathr_edges(df_dict, perm=-1, p_thresh_pos=None, p_thresh_neg=None, r_thresh_pos=None, r_thresh_neg=None, percentile_neg=None, percentile_pos=None, top_n_pos=None, top_n_neg=None):
  folds_list = list(df_dict[perm].keys())
  n_folds = len(folds_list)
  n_edges = len(df_dict[perm][folds_list[0]])
  all_masks = {}
  all_masks['pos'] = np.zeros((n_folds, n_edges))
  all_masks['neg'] = np.zeros((n_folds, n_edges))
  
  for fold in folds_list:
    pcorr_df = df_dict[perm][fold]
    suprathr_edges_mask = {}
    if p_thresh_pos and p_thresh_neg:
      suprathr_edges_mask['pos'] = (pcorr_df['r'] > 0) & (pcorr_df['p-val'] <= p_thresh_pos)
      suprathr_edges_mask['neg'] = (pcorr_df['r'] < 0) & (pcorr_df['p-val'] <= p_thresh_neg)
    elif r_thresh_pos and r_thresh_neg:
      suprathr_edges_mask['pos'] = pcorr_df['r'] > r_thresh_pos
      suprathr_edges_mask['neg'] = pcorr_df['r'] < -abs(r_thresh_neg) # r_thresh_neg can be both given as a positive or a negative value
    elif percentile_pos and percentile_neg:
      r_thresh_pos = np.nanpercentile(pcorr_df['r'], percentile_pos)
      r_thresh_neg = np.nanpercentile(pcorr_df['r'][pcorr_df['r'] < 0], 100 - percentile_neg)
      suprathr_edges_mask['pos'] = pcorr_df['r'] > r_thresh_pos
      suprathr_edges_mask['neg'] = pcorr_df['r'] < -abs(r_thresh_neg)
    elif top_n_pos and top_n_neg:
      suprathr_edges_mask['pos'] = np.zeros(pcorr_df.shape[0])
      suprathr_edges_mask['neg'] = np.zeros(pcorr_df.shape[0])
      suprathr_edges_mask['pos'][np.argpartition(pcorr_df['r'][pcorr_df['r'].notna()], -top_n_pos)[-top_n_pos:]] = 1
      suprathr_edges_mask['neg'][np.argpartition(pcorr_df['r'][pcorr_df['r'].notna()], top_n_neg)[:top_n_neg]] = 1
    else:
      raise TypeError('Either p_thresh_{neg, pos} or r_thresh_{neg, pos} or percentile_{neg, pos} or top_n_{pos, neg} needs to be defined.')
    
    printv("Fold {}: Pos/neg suprathreshold edges (max r pos/max r neg): {}/{} ({}/{})".format(fold+1, suprathr_edges_mask['pos'].sum(), suprathr_edges_mask['neg'].sum(), pcorr_df['r'].max(), pcorr_df['r'].min()))
    all_masks['pos'][fold,:] = suprathr_edges_mask['pos'].astype(bool)
    all_masks['neg'][fold,:] = suprathr_edges_mask['neg'].astype(bool)
  
  return all_masks

def start_autosave(path, ray_handler, save_size=1000):
    from ray.util.ml_utils.node import force_on_current_node
    AutoSaveActor = force_on_current_node(AutoSaveActor)
    autosave_actor = AutoSaveActor.remote(path, ray_handler, save_size)

    return autosave_actor

class RayHandler:
  def __init__(self, fc_data, behav_data, behav, covars, n_perm=0, **ray_kwargs):
    self.behav_data = behav_data # For adding kfold_indices
    self.fselection_results = {}
    self.fselection_results[-1] = {} # Create sub-dictionary for original (i.e. non-permuted) data
    self.prediction_results = {}

    ray.shutdown() # Make sure Ray is only initialised once
    self.ray_info = ray.init(**ray_kwargs)

    self.actors_list = []
    self.job_list = []

    self.data_dict = {}
    self.data_dict['behav'] = behav
    self.data_dict['covars'] = covars
    self.data_dict['n_perm'] = n_perm
    self.data_dict['data'] = fc_data
    self.data_dict['edges'] = self.data_dict['data'].columns.astype(str) # Save edges columns before adding behavioral columns
    # Pingouin needs all the data (edges, behav, and covars) in a single DataFrame
    if covars:
      self.data_dict['data'][covars] = behav_data[covars]
    if n_perm > 0:
      # It seems to be more efficient to create a separate df and concat later;
      # .to_frame() converts Pandas series into a DataFrame on-the-fly
      behav_df = behav_data[behav].to_frame()
      for perm in range(n_perm):
        behav_df["{}-perm-{}".format(behav, perm)] = np.random.permutation(behav_df[behav])
      behav_df = behav_df.copy()
      # To avaid fragmentation (and the corresponding warning), consolidate into a 
      # new DataFrame)
      self.data_dict['data'] = pd.concat([self.data_dict['data'], behav_df], axis=1)
    else:
      self.data_dict['data'][behav] = behav_data[behav]
    self.data_dict['data'].columns = self.data_dict['data'].columns.astype(str)

    # Start status and results actor preferrably (status actor) or forcibly
    # (results actor) on head node
    self.status_actor = StatusActor.options(
      scheduling_strategy=ray.util.scheduling_strategies.NodeAffinitySchedulingStrategy(
      node_id = ray.get_runtime_context().node_id,
      soft = True)).remote()
    self.results_actor = ResultsActor.options(
      scheduling_strategy=ray.util.scheduling_strategies.NodeAffinitySchedulingStrategy(
      node_id = ray.get_runtime_context().node_id,
      soft = False)).remote()
  
  def add_kfold_indices(self, n_folds, clean=True):
    subject_ids = self.data_dict['data'].index
    kfold_indices = get_kfold_indices(subject_ids, n_folds)
    if clean:
        kfold_indices = clean_kfold_indices(kfold_indices, self.behav_data)
    self.data_dict['kfold_indices'] = kfold_indices
    printv("You need to (re-) upload data after this operation.")
    
  def upload_data(self):
    self.data_object = ray.put(self.data_dict)

  def start_actors(self, job_list):
    printv("Starting actors for {} jobs ...".format(len(job_list)))
    self.job_list.extend(job_list)
    self.actors = [RayActor.remote(self.data_object, job, self.status_actor, self.results_actor) for job in job_list]
    
  def get_fselection_results(self):
      results = ray.get(self.results_actor.get_fselection_results.remote())
      n = 1
      N = len(results)
      for result in results:
          fold = result[0]
          perm = result[1]
          df = result[2]
          self.fselection_results[perm][fold] = df
          n += 1
      return self.fselection_results

  def get_prediction_results(self):
      results = ray.get(self.results_actor.get_prediction_results.remote())
      for results_dict in results:
          if results_dict['perm'] not in self.prediction_results:
              self.prediction_results[results_dict['perm']] = pd.DataFrame()
              self.prediction_results[results_dict['perm']]['observed'] = self.data_dict['data'][self.data_dict['behav']]
          for tail in ('pos', 'neg', 'glm'):
              self.prediction_results[results_dict['perm']].loc[results_dict['test_IDs'], [tail]] = results_dict[tail]
      return self.prediction_results

  def status(self, verbose=True):
    n = 1
    self.status_dict = ray.get(self.status_actor.get_status.remote())
    for pid, info in self.status_dict.items():
        if(info['msg']): # Only print alive actors (-> msg != None)
            print("Actor {} [{}@{}]: {}".format(n, pid, info['node'], info['msg']))
            n += 1
    print("\n")
    n_done = ray.get(self.results_actor.get_size.remote())
    n_remaining = len(self.job_list) - n_done
    print("Jobs done: {}".format(n_done))
    print("Jobs remaining: {}".format(n_remaining))
    
    return n_done, n_remaining
        
  def terminate(self):
    ray.shutdown()

@ray.remote
class StatusActor:
    def __init__(self):
        self.status_dict = {}

    def send(self, status):
        pid = status[0]
        node = status[1]
        msg = status[2]
        self.status_dict[pid] = {"msg": msg, "node": node}

    def exit(self):
        ray.actor.exit_actor()

    def get_status(self):
        return self.status_dict

@ray.remote
class ResultsActor:
    def __init__(self):
        # Create dictionaries to keep results (it makes sense to do this class-wide to add results on-the-fly
        # and for later reference if get results functions are called too early for example
        self.fselection_results = []
        self.prediction_results = []

    def send(self, results):
        """
        Receives result from worker, determines type of result and calls the appropriate
        processing function
        """
        if type(results) == list:
            self.fselection_results.append(results)
        elif type(results) == dict:
            self.prediction_results.append(results)
        else: # Make sure no results are lost
            self.fselection_results.append(results)

    def get_fselection_results(self):
        return self.fselection_results

    def get_prediction_results(self):
        return self.prediction_results

    def save_prediction_results(self, path):
        np.save(path, self.prediction_results)

    def get_size(self):
        size = len(self.fselection_results) + len(self.prediction_results)
        return size

    def exit(self):
        ray.actor.exit_actor()

@ray.remote(num_cpus=1)
class RayActor:
    def __init__(self, data, job, status_actor, results_actor):
        self.data = data
        self.status_actor = status_actor
        self.results_actor = results_actor

        self.node = socket.gethostname()
        self.pid = os.getpid()
    
        self.start_job(job)

    def exit(self):
        self.status_update(None)
        ray.actor.exit_actor()
    
    def start_job(self, job):
        job_type = job[0]
        fold = job[1]
        perm = job[2]
        try:
            obj = job[3] # This was added to optionally supply train_subs manually (for debugging or learning purposes); when predicting, a mask _has_ to be supplied
        except IndexError:
            obj = None
    
        if job_type == 'fselection':
            fselection_result = self.do_fselection(fold, perm, train_subs=obj)
            self.results_actor.send.remote([fold, perm, fselection_result])
            self.exit() # Exit so memory gets freed up and no substantial memory leak happens
        elif job_type == 'prediction':
            prediction_result = self.do_prediction(fold, perm, mask=obj)
            self.results_actor.send.remote(prediction_result)
            self.exit()
        elif job_type == 'fselection_and_prediction':
            fselection_result = self.do_fselection(fold, perm, train_subs=obj)
            df_dict = {fold: fselection_result}
            self.status_update("Thresholding edges...")
            mask = self.get_suprathr_edges_new(df_dict, p_thresh_pos=0.001, p_thresh_neg=0.001) # TODO: Make this configurable via **kwargs or something (or as an option in data_dict?)
            prediction_result = self.do_prediction(fold, perm, mask=mask[fold])
            self.results_actor.send.remote(prediction_result)
            self.exit()
        else:
            raise TypeError('Ill-defined job type.')
        
    def get_suprathr_edges_new(self, df_dict, p_thresh_pos=None, p_thresh_neg=None, r_thresh_pos=None, r_thresh_neg=None, percentile_neg=None, percentile_pos=None, top_n_pos=None, top_n_neg=None):
        # We need to duplicate this here or "no module found" error will occur
        folds_list = list(df_dict.keys())
        n_edges = len(df_dict[folds_list[0]])
        masks_dict = {}

        for fold in folds_list:
            pcorr_df = df_dict[fold]
            n_edges = len(df_dict[fold])
            masks_dict[fold] = {}
            suprathr_edges_mask = {}
            if p_thresh_pos and p_thresh_neg:
                suprathr_edges_mask['pos'] = (pcorr_df['r'] > 0) & (pcorr_df['p-val'] <= p_thresh_pos)
                suprathr_edges_mask['neg'] = (pcorr_df['r'] < 0) & (pcorr_df['p-val'] <= p_thresh_neg)
            elif r_thresh_pos and r_thresh_neg:
                suprathr_edges_mask['pos'] = pcorr_df['r'] > r_thresh_pos
                suprathr_edges_mask['neg'] = pcorr_df['r'] < -abs(r_thresh_neg) # r_thresh_neg can be both given as a positive or a negative value
            elif percentile_pos and percentile_neg:
                r_thresh_pos = np.nanpercentile(pcorr_df['r'], percentile_pos)
                r_thresh_neg = np.nanpercentile(pcorr_df['r'][pcorr_df['r'] < 0], 100 - percentile_neg)
                suprathr_edges_mask['pos'] = pcorr_df['r'] > r_thresh_pos
                suprathr_edges_mask['neg'] = pcorr_df['r'] < -abs(r_thresh_neg)
            elif top_n_pos and top_n_neg:
                suprathr_edges_mask['pos'] = np.zeros(pcorr_df.shape[0])
                suprathr_edges_mask['neg'] = np.zeros(pcorr_df.shape[0])
                suprathr_edges_mask['pos'][np.argpartition(pcorr_df['r'][pcorr_df['r'].notna()], -top_n_pos)[-top_n_pos:]] = 1
                suprathr_edges_mask['neg'][np.argpartition(pcorr_df['r'][pcorr_df['r'].notna()], top_n_neg)[:top_n_neg]] = 1
            else:
                raise TypeError('Either p_thresh_{neg, pos} or r_thresh_{neg, pos} or percentile_{neg, pos} or top_n_{pos, neg} needs to be defined.')

            for tail in ('pos', 'neg'):
                masks_dict[fold][tail] = np.zeros(n_edges)
                masks_dict[fold][tail][:] = suprathr_edges_mask[tail].astype(bool)

        return masks_dict
  

    def do_fselection(self, fold, perm, train_subs=None):
        if not train_subs: # Get train_subs from uploaded database if not supplied
            train_subs = self.data['kfold_indices']['train'][fold]
        fselection_result = self.edgewise_pcorr(train_subs, fold, perm)
        
        return fselection_result
  
    def do_prediction(self, fold, perm, mask=None):
        if not mask:
            raise ValueError("A mask needs to be supplied in the job description.")
        train_subs = self.data['kfold_indices']['train'][fold]
        test_subs = self.data['kfold_indices']['test'][fold]
        self.status_update("Predicting fold {}...".format(fold+1)) # No detailed update needed, so once is sufficient
        prediction_result = self.predict(mask, train_subs, test_subs, perm)
      
        return prediction_result
    
    def status_update(self, msg):
        self.status_actor.send.remote([self.pid, self.node, msg])
    
    def edgewise_pcorr(self, train_subs, fold, perm, method='pearson'):
        corr_dfs = [] # Appending to list and then creating dataframe is substantially faster than appending to dataframe
        empty_df = pd.DataFrame({'r': {'pearson': np.nan}, 'p-val': {'pearson': np.nan}}) # Handle all-zero edges
        train_data = self.data['data'].loc[train_subs]
        train_data.columns = train_data.columns.astype(str)
    
        N = len(self.data['edges'])
        n = 1
        percent = round((n/N)*100)
    
        for edge in self.data['edges']: # Edge-wise correlation
            if (train_data[edge] != 0).any(): # All-zero columns will raise a ValueError exception. This is _way_ faster than try: except:
                if perm >= 0:
                    y = "{}-perm-{}".format(self.data['behav'], perm)
                else:
                    y = self.data['behav']
                if self.data['covars']:
                    pcorr = partial_corr(data=train_data, x=edge, y=y, covar=self.data['covars'], method=method)[['r', 'p-val']] # Taking only the necessary columns speeds this up a few %
                    pcorr['covars'] = True # Debug, remove later
                else: # We could also use pcorr from Pingouin on the entire df, but this is a prohibitively memory-intensive operation; edge-wise like this works just fine. This was introduced to test implausibly good results for the above operation by Pengouin's partial_corr, by setting covars=None we can use SciPy's implementation of Pearson's r
                    pcorr = empty_df.copy()
                    pcorr[['r', 'p-val']] = sp.stats.pearsonr(train_data.loc[:, edge], train_data.loc[:, y]) # We are basically reproducing Pingouin's output format here for unified downstream processing
                    pcorr['covars'] = False # Debug, remove later
            else:
                pcorr = empty_df
            corr_dfs.append(pcorr)
            percent_new = round((n/N)*100)
            if perm >= 0:
                fold_msg = "{} of permutation {}".format(fold+1, perm+1)
            else:
                fold_msg = fold+1
            if percent_new > percent:
                self.status_update("Computing fold {} ({} %)...".format(fold_msg, percent_new))
                percent = percent_new
            n += 1
        self.status_update("Assembling data frame...")
        combined_corr_dfs = pd.concat(corr_dfs) # Assembling df before .put() seems to avoid awfully slow pickling of data through queue (or whatever, it is orders of magnitude faster that way)
        return combined_corr_dfs

    
    def predict(self, mask_dict, kfold_indices_train, kfold_indices_test, perm):
        train_vcts = pd.DataFrame(self.data['data'].loc[kfold_indices_train, self.data['edges']])
        test_vcts = pd.DataFrame(self.data['data'].loc[kfold_indices_test, self.data['edges']])
        
        if perm >= 0:
            behav = "{}-perm-{}".format(self.data['behav'], perm)
        else:
            behav = self.data['behav']
            
        train_behav = pd.DataFrame(self.data['data'].loc[kfold_indices_train, behav])
    
        model = self.build_models(mask_dict, train_vcts, train_behav)
        behav_pred = self.apply_models(mask_dict, test_vcts, model)
        behav_pred["IDs"] = kfold_indices_test
        behav_pred["model"] = model # Debug
        behav_pred["mask"] = mask_dict # Debug
        behav_pred["train_IDs"] = kfold_indices_train # Debug
        behav_pred["test_IDs"] = kfold_indices_test # Debug
        behav_pred["perm"] = perm # Debug
    
        return behav_pred
        
    
    def build_models(self, mask_dict, train_vcts, train_behav):
        """
        Takes a feature mask, sums all edges in the mask for each subject, and uses simple linear 
        regression to relate summed network strength to behavior; returns a dictionary with the    
        model
        """

        assert train_vcts.index.equals(train_behav.index), "Row indices of FC vcts and behavior don't match!"
        lr = LinearRegression()
        model_dict = {}

        X_glm = np.zeros((len(train_vcts.index), len(mask_dict.items())))

        t = 0
        # Loop through pos and neg tails
        for tail, mask in mask_dict.items():
            summed_edges = train_vcts.values[:, mask.astype(bool)].sum(axis=1)
            X_glm[:, t] = summed_edges
            lr.fit(summed_edges.reshape(-1, 1), train_behav)
            model_dict[tail] = (lr.coef_, lr.intercept_)
            t += 1
        
        X_glm = np.c_[X_glm, np.ones(X_glm.shape[0])]
        model_dict["glm"] = tuple(np.linalg.lstsq(X_glm, train_behav, rcond=None)[0])

        return model_dict

    def apply_models(self, mask_dict, test_vcts, model_dict):
        """
        Applies a previously trained linear regression model to a test set to generate 
        predictions of behavior.
        """

        behav_pred = {}

        X_glm = np.zeros((test_vcts.shape[0], len(mask_dict.items())))

        # Loop through pos and neg tails
        t = 0
        for tail, mask in mask_dict.items():
            X = test_vcts.loc[:, mask.astype(bool)].sum(axis=1)
            X_glm[:, t] = X

            slope, intercept = model_dict[tail]
            behav_pred[tail] = slope[0]*X + intercept
            t+=1

        X_glm = np.c_[X_glm, np.ones(X_glm.shape[0])]
        glm_results = np.dot(X_glm, model_dict["glm"])
        behav_pred["glm"] = glm_results[:, 0] # Transforms 2d array into 1d array

        return behav_pred
