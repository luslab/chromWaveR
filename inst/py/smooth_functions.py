import numpy as np
from copy import copy
from scipy.ndimage import gaussian_filter1d

def int_to_float(x, u=1):
  # scale if necessary
  x =x / np.float(u)
  return x

### ChromWAVE model smoothing functions
# y are the predictions of the model. eg
# load model
# chromwave_model = enomic_wavenet.GenomeWaveNet()
# chromwave_model.deserialize(directory = TF_Nuc_model_dir)
# first change input lenght to be able to predict more than 4000bp
# chromwave_model.change_input_length()
# x are the 1-hot encoded chromosomes
# y_p = chromwave_model.predict(x)
# y= [y.argmax(axis=2) for y in y_p]
def invert_discretizing(y_p,preprocessing_parameters):
  y = [y.argmax(axis=2) for y in y_p]
  y_smooth = [_smooth_discrete_data(copy(seq), preprocessing_params=params)
              for (seq, params) in zip(y, preprocessing_parameters)]
  return [y-np.mean(y) for y in y_smooth]

def _smooth_discrete_data(y,preprocessing_params):
  discretize_function=preprocessing_params['discretize_function']
  if discretize_function is not None:
    if discretize_function == 'float_to_int':
      u = preprocessing_params['u']
      y = [np.array(int_to_float(seq, u=u)) for seq in y]
      y = np.vstack([np.array(gaussian_filter1d(seq, sigma=3,truncate=4.0)) for seq in y])
      return y
    else:
      print('Smoothing not implememented yet.')

def get_preprocessing_parameters(model_type = 'TF-Nuc'):
  if model_type == 'TF-Nuc':
    # TF signal smoothing parameters
    tf_preprocessing_params = {'times_median_coverage_max': 3,
                               'discretize_function': 'float_to_int',
                               'assign_non_covered_bases': None,
                               'u': 10, 'smooth_signal': True,
                               'discretize_function': 'float_to_int',
                               'sigma': 5, 'truncate': 3,
                               'smoothing_function': 'gaussian_filter1d',
                               'x_thresh': 0, 'run_thresh': 10,
                               'normalise_read_counts': 'genome_mean'}
    # Nuc signal smoothing parameters
    nuc_preprocessing_params = {'times_median_coverage_max': 3,
                                'discretize_function': 'float_to_int',
                                'assign_non_covered_bases': None,
                                'u': 5, 'smooth_signal': True,
                                'discretize_function': 'float_to_int',
                                'sigma': 5, 'truncate': 3,
                                'smoothing_function': 'gaussian_filter1d',
                                'x_thresh': 0, 'run_thresh': 50,
                                'normalise_read_counts': 'genome_mean'}
    return([tf_preprocessing_params, nuc_preprocessing_params])
  elif (model_type == 'invivo' or model_type == 'invitro'):
    # Yeast invivo/invitro smoothing parameters
    preprocessing_params = {'times_median_coverage_max': 3,
                            'discretize_function': 'float_to_int',
                            'assign_non_covered_bases': None,
                            'u': 2.5, 'smooth_signal': True,
                            'sigma': 5, 'truncate': 3,
                            'smoothing_function': 'gaussian_filter1d',
                            'x_thresh': 0, 'run_thresh': 50,
                            'normalise_read_counts': None}
    return([preprocessing_params])
  elif model_type == 'hs_promoter':
    # Human promoter nuc smoothing parameters
    preprocessing_params = {'times_median_coverage_max': 3,
                            'discretize_function': 'float_to_int',
                            'assign_non_covered_bases': 'chrom_mean',
                            'u': 2, 'smooth_signal': True,
                            'discretize_function': 'float_to_int',
                            'sigma': 5, 'truncate': 3,
                            'smoothing_function': 'gaussian_filter1d',
                            'x_thresh': 0, 'run_thresh': 50,
                            'normalise_read_counts': 'log2_ratio'}
    return([preprocessing_params])



