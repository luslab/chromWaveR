import numpy as np
import re
from keras.models import model_from_json

def load_chromwave_model(json_file):
  # Load model from *.json file and change input shape
  json_file = open(json_file, 'r')
  model = json_file.read()
  old_input = '\"batch_input_shape\": \[null, [0-9]{1,9}, 4\]'
  new_input = '"batch_input_shape": [null, null, 4]'
  model = re.sub(old_input, new_input, model)
  json_file.close()
  # Compile model
  model = model_from_json(model)
  return model

def one_hot_encode(sequences):
  seqs = [np.transpose(bases_to_binary(seq)) for seq in sequences]
  return seqs

def bases_to_binary(sequence, trim_start=0, trim_end=0, paddingL=0,paddingR=0 , matrix_bases='ACGT'):

    binarised_bases = np.array(
        [
            [# trim start and end characters if necessary and force upper case for sequence
                (x == base) for x in list(sequence[trim_start:(None if trim_end == 0 else -trim_end)].upper())  # force upper case for bases too
            ] for base in matrix_bases.upper()
        ],dtype=np.int8
    )

    # if needed, add padding of Ns on either side of sequence and return
    if paddingL > 0 or paddingR>0:
        return padding_dna(binarised_bases, paddingL, paddingR)
    # otherwise just return the binarised bases
    return binarised_bases
