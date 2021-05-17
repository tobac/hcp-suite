import sys
sys.path.append("/home/sc.uni-leipzig.de/ty046wyna/hcp-suite/lib")
import ray
from ray.util.queue import Queue
from cpm_cov import *
import numpy as np
import pandas as pd
from time import sleep
import os
import pickle
from sklearn.model_selection import KFold
from pingouin import partial_corr
from sklearn.linear_model import LinearRegression
from nilearn import plotting
from hcpsuite import load_list_from_file, get_nimg_data, printv, nifti_dim_to_cifti_dim, symmetrize_matrices, symmetrize, timer

def get_behav_data(fname, ids):
  """
  Read behavioural data for specified IDs from CSV and return dataframe
  """
  
  behav_data = pd.read_csv(fname, dtype={'Subject': str}, low_memory = False)
  behav_data.set_index('Subject', inplace=True)
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

def plot_consistent_edges(all_masks, tail, thresh = 1., color='gray', coords=None):
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
                    edge_kwargs={"linewidth": 1, 'color': color})


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

def get_suprathr_edges(df_dict, p_thresh_pos=None, p_thresh_neg=None, r_thresh_pos=None, r_thresh_neg=None, percentile_neg=None, percentile_pos=None):
  n_folds = len(df_dict)
  n_edges = len(df_dict[0])
  all_masks = {}
  all_masks['pos'] = np.zeros((n_folds, n_edges))
  all_masks['neg'] = np.zeros((n_folds, n_edges))
  
  for fold in range(n_folds):
    pcorr_df = df_dict[fold]
    suprathr_edges_mask = {}
    if p_thresh_pos and p_thresh_neg:
      suprathr_edges_mask['pos'] = (pcorr_df['r'] > 0) & (pcorr_df['p-val'] <= p_thresh_pos)
      suprathr_edges_mask['neg'] = (pcorr_df['r'] < 0) & (pcorr_df['p-val'] <= p_thresh_neg)
    elif r_thresh_pos and r_thresh_neg:
      suprathr_edges_mask['pos'] = pcorr_df['r'] > r_thresh_pos
      suprathr_edges_mask['neg'] = pcorr_df['r'] < -abs(r_thresh_neg) # r_thresh_neg can be both given as a positive or a negative value
    elif percentile_pos and percentile_neg:
      r_thresh_pos = np.percentile(pcorr_df['r'], percentile_pos)
      r_thresh_neg = np.percentile(pcorr_df['r'][pcorr_df['r'] < 0], 100 - percentile_neg)
      suprathr_edges_mask['pos'] = pcorr_df['r'] > r_thresh_pos
      suprathr_edges_mask['neg'] = pcorr_df['r'] < -abs(r_thresh_neg)  
    else:
      raise TypeError('Either p_thresh_{neg, pos} or r_thresh_{neg, pos} or percentile_{neg, pos} needs to be defined.')
    
    printv("Fold {}: Pos/neg suprathreshold edges (max r pos/max r neg): {}/{} ({}/{})".format(fold+1, suprathr_edges_mask['pos'].sum(), suprathr_edges_mask['neg'].sum(), pcorr_df['r'].max(), pcorr_df['r'].min()))
    all_masks['pos'][fold,:] = suprathr_edges_mask['pos'].astype(bool)
    all_masks['neg'][fold,:] = suprathr_edges_mask['neg'].astype(bool)
  
  return all_masks

class RayHandler:
  def __init__(self, fc_data, behav_data, behav, covars, ray_address=None):
    ray.shutdown() # Make sure Ray is only initialised once
    if ray_address:
      self.ray_info = ray.init(address='auto')
    else:
      self.ray_info = ray.init()

    self.in_queue = Queue()
    self.out_queue = Queue()
    self.status_queue = Queue()
    
    self.status_dict = {}

    self.data_dict = {}
    self.data_dict['behav'] = behav
    self.data_dict['covars'] = covars
    self.data_dict['data'] = fc_data
    self.data_dict['edges'] = self.data_dict['data'].columns.astype(str) # Save edges columns before adding behavioral columns
    if covars:
      self.data_dict['data'][covars] = behav_data[covars]
    self.data_dict['data'][behav] = behav_data[behav]
    self.data_dict['data'].columns = self.data_dict['data'].columns.astype(str)
    
  def upload_data(self):
    # Allows us to manipulate data in-class before uploading
    # TODO: Put this and start_workers function in __init__() again?
    self.data_object = ray.put(self.data_dict)
    
  def start_workers(self, n_workers):
    printv("Starting {} workers".format(n_workers))
    self.workers = [RayWorker.remote(self.data_object, self.in_queue, self.out_queue, self.status_queue) for _ in range(n_workers)]
    
  def submit_fselection(self, train_subs, fold):
    self.in_queue.put(['fselection', train_subs, fold])
  
  def submit_prediction(self, mask, kfold_indices_train, kfold_indices_test, fold):
    self.in_queue.put(['prediction', mask, kfold_indices_train, kfold_indices_test, fold])
    
  def get_fselection_results(self):
    results = {}
    n = 1
    N = self.out_queue.qsize()
    while not self.out_queue.empty():
        printv("Retrieving item {} of {}".format(n, N), update=True)
        result = self.out_queue.get()
        results[result[0]] = result[1]
        n += 1
    return results

  def get_prediction_results(self):
    results = pd.DataFrame()
    results['observed'] = self.data_dict['data'][self.data_dict['behav']]
    for tail in ('pos', 'neg', 'glm'): # Initialize rows
        results[tail] = np.nan
    
    n = 1
    N = self.out_queue.qsize()
    while not self.out_queue.empty():
        printv("Retrieving item {} of {}".format(n, N), update=True)
        result_dict = self.out_queue.get()
        for tail in ('pos', 'neg', 'glm'):
          results.loc[result_dict['IDs'], tail] = result_dict[tail]
        n += 1
    return results

  def status(self):
    n = 1
    N = self.status_queue.qsize()
    while not self.status_queue.empty():
      printv("Retrieving item {} of {}".format(n, N), update=True)
      status_list = self.status_queue.get()
      pid = status_list[0]
      msg = status_list[1]
      self.status_dict[pid] = msg
      n += 1
    n = 1
    N_workers = len(self.status_dict)
    width = len(str(N_workers))
    for pid, msg in self.status_dict.items():
        print("Worker {}/{} [{}]: {}".format(str(n).zfill(width), N_workers, pid, msg))
        n += 1
    print("\n")
    print("Folds remaining in queue: {}".format(self.in_queue.qsize()))
        
  def terminate(self):
    ray.shutdown()

@ray.remote
class RayWorker:
  def __init__(self, data, in_queue, out_queue, status_queue):
    
    self.data = data
    self.behav = self.data['behav']
    self.covars = self.data['covars']
    
    self.in_queue = in_queue
    self.out_queue = out_queue
    self.status_queue = status_queue

    self.pid = os.getpid()
    
    self.in_queue_listener() # Listen for jobs
    
  def in_queue_listener(self):
    while True:
      self.status_update("Listening for jobs...")
      while self.in_queue.empty():
        sleep(0.5)
      self.status_update("Receiving job from queue...")
      obj_list = self.in_queue.get()
      job_type = obj_list[0]
    
      if job_type == 'fselection':
        train_subs = obj_list[1]
        fold = obj_list[2]
        self.edgewise_pcorr(train_subs, fold) 
      elif job_type == 'prediction':
        mask = obj_list[1]
        kfold_indices_train = obj_list[2]
        kfold_indices_test = obj_list[3]
        fold = obj_list[4]
        self.status_update("Predicting fold {}...".format(fold+1)) # No detailed update needed, so once is sufficient
        self.predict(mask, kfold_indices_train, kfold_indices_test)
      else:
        raise TypeError('Ill-defined job type.')
  
  def status_update(self, msg):
    self.status_queue.put([self.pid, msg])
        
  def edgewise_pcorr(self, train_subs, fold, method='pearson'):
    corr_dfs = [] # Appending to list and then creating dataframe is substantially faster than appending to dataframe
    empty_df = pd.DataFrame({'r': {'pearson': np.nan}, 'p-val': {'pearson': np.nan}}) # Handle all-zero edges
    train_data = self.data['data'].loc[train_subs]
    train_data.columns = train_data.columns.astype(str)
    
    N = len(self.data['edges'])
    n = 1
    percent = round((n/N)*100)
    
    for edge in self.data['edges']: # Edge-wise correlation
      if (train_data[edge] != 0).any(): # All-zero columns will raise a ValueError exception. This is _way_ faster than try: except:
        if self.covars:
          pcorr = partial_corr(data=train_data, x=edge, y=self.behav, covar=self.covars, method=method)[['r', 'p-val']] # Taking only the necessary columns speeds this up a few %
          pcorr['covars'] = True # Debug, remove later
        else: # We could also use pcorr from Pengouin on the entire df, but this is a prohibitively memory-intensive operation; edge-wise like this works just fine. This was introduced to test implausibliy good results for the above operation by Pengouin's partial_corr, by setting covars=None we can use SciPy's implementation of Pearson's r
          pcorr = empty_df.copy()
          pcorr[['r', 'p-val']] = sp.stats.pearsonr(train_data.loc[:, edge], train_data.loc[:, self.behav]) # We are basically reproducing Pengouin's output format here for unified downstream processing
          pcorr['covars'] = False # Debug, remove later
      else:
        pcorr = empty_df
      corr_dfs.append(pcorr)
      percent_new = round((n/N)*100)
      if percent_new > percent:
        self.status_update("Computing fold {} ({} %)...".format(fold+1, percent_new))
        percent = percent_new
      n += 1
    self.status_update("Assembling data frame...")
    self.out_queue.put([fold, pd.concat(corr_dfs)]) # Assembling df before .put() seems to avoid awfully slow pickling of data through queue (or whatever, it is orders of magnitude faster that way)
    self.status_update("Listening for jobs...")

  def predict(self, mask_dict, kfold_indices_train, kfold_indices_test):
    train_vcts = pd.DataFrame(self.data['data'].loc[kfold_indices_train, self.data['edges']])
    test_vcts = pd.DataFrame(self.data['data'].loc[kfold_indices_test, self.data['edges']])
    train_behav = pd.DataFrame(self.data['data'].loc[kfold_indices_train, self.data['behav']])
    
    model = self.build_models(mask_dict, train_vcts, train_behav)
    behav_pred = self.apply_models(mask_dict, test_vcts, model)
    behav_pred["IDs"] = kfold_indices_test
    behav_pred["model"] = model # Debug
    behav_pred["mask"] = mask_dict # Debug
    behav_pred["train_IDs"] = kfold_indices_train # Debug
    behav_pred["test_IDs"] = kfold_indices_test # Debug
    self.out_queue.put(behav_pred)
        
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
    behav_pred["glm"] = np.dot(X_glm, model_dict["glm"])

    return behav_pred