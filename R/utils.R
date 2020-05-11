#' List all available chromWave models.
#'
#' \code{availableModels} returns all previously trained and available
#' chromWave nucleosome models.
#'
#' \tabular{lll}{
#'  \strong{Model} \tab \strong{Organism} \tab \strong{Description}\cr
#'   invitro \tab sacCer \tab Model was trained on in vitro nucleosome data from Kaplan et al (2009) \cr
#'   invivo \tab sacCer \tab Model was trained on in vivo nucleosome data from Kaplan et al (2009) \cr
#'   TF-NUC \tab sacCer \tab Model was trained on MNase-seq data from Henikoff et al (2011)\cr
#'   hs_promoter \tab hs \tab Model was trained on promoter subset of MNase-seq data from Gaffney et al (2012) \cr
#' }
#'
#' @import tibble
#'
#' @return Function returns tibble of all chromwaveR included models.
#'
#' @export
#'
#' @examples
#' availableModels()
#'
availableModels <- function() {
  tibble(model_name = c('invitro', 'invivo', 'TF-Nuc', 'hs_promoter'),
         organism = c('sacCer', 'sacCer', 'sacCer', 'hs'),
         description = c('Model was trained on in vitro nucleosome data from Kaplan et al (2009)',
                         'Model was trained on in vivo nucleosome data from Kaplan et al (2009)',
                         'Model was trained on MNase-seq data from Henikoff et al (2011)',
                         'Model was trained on promoter subset of MNase-seq data from Gaffney et al (2012)'))
}


#' Load chromWave nucleosome model.
#'
#' \code{loadChromWaveModel} loads pre-trained nucleosome model from *.json/*.h5
#' files included in this package. The model input shape is altered by this function
#' to (None, None, 4) enabling chromosome wide predicitons.
#' Function serves as wrapper for python function load_chromwave_model.
#' Please check \code{availableModels()} for list of pretrained nucleosome models.
#'
#' @param \code{model.name} Models include invitro, invivo, TF-NUC and hs_promoter.
#'
#' @import tibble
#' @importFrom reticulate source_python
#' @importFrom kerasR keras_load_weights
#'
#' @return Function returns chromWave model as keras object.
#'
#' @export
#'
#' @examples
#' ## Load yeast in vivo nucleosome model
#' model <- loadChromWaveModel(model_name = 'invivo')
#'
#' ## Load human promoter nucleosome model
#' model <- loadChromWaveModel(model_name = 'hs_promoter')
#'
loadChromWaveModel <- function(model.name) {
  # Check for valid model.name
  if(!(model.name %in% availableModels()$model_name)) {
    stop('Invalid model.name.
         Please check availableModels() for valid model names.')
  }

  # Model path
  model.path <- system.file('models', package = 'chromWaveR')

  # Find model associated *.json file
  json.file <- list.files(model.path, pattern = sprintf('%s.*.json', model.name),
                          full.names = T)
  h5.file <- list.files(model.path, pattern = sprintf('%s.*.h5', model.name),
                        full.names = T)

  # Load python funcitons: load_chromwave_model
  utils.py <- system.file('py', 'utils.py', package = 'chromWaveR')
  source_python(utils.py)

  # Load chromwave nucleosome model
  model <- load_chromwave_model(json.file)
  keras_load_weights(model, h5.file)

  return(model)
}

#' Function to create random DNA strings.
#'
#' \code{rmdDNA} creates a DNAStringSet with \code{n} DNAStrings with
#' the width \code{width}.
#'
#' @param \code{n} Number of rmd generated DNA strings.
#' @param \code{width} Length in bp of generated DNA strings.
#'
#' @importFrom Biostrings DNAString DNAStringSet
#'
#' @return Function returns randomly generated DNAStringSet.
#'
#' @export
#'
#' @examples
#' # Create 100 rmd DNA sequences of 1000bp lenght
#' rmdDNA(n = 100, width = 10^3)
#'
rmdDNA <- function(n = 10, width = 500) {
  # Sample nucleotides and create DNA strings
  dna.strs <- lapply(1:n, function(i) {
    dna.str <- sample(c('A', 'C', 'G', 'T'), size = width, replace = T)
    dna.str <- DNAString(paste(dna.str, collapse = ''))
    return(dna.str)
  })
  dna.strs <- DNAStringSet(dna.strs)
  return(dna.strs)
}

#' Convert DNA sequence from string to one hot encoded matrix.
#'
#' \code{oneHotEncode} is a wrapper function for the
#' chromWave python function one_hot_encode.
#' Function converts DNA strings to one-hot encoded sequence matrices.
#'
#' @param \code{dna.strs} DNAStringSet
#'
#' @importFrom reticulate source_python
#'
#' @return Function returns list of one-hot encoded DNAStrings.
#'
#' @export
#'
#' @examples
#' # Convert DNAStringSet to list of one-hot encoded matrices
#' dna.strs <- rmdDNA()
#' oneHot.strs <- oneHotEncode(dna.strs)
#'
oneHotEncode <- function(dna.strs) {
  # Load additional python funcitons
  utils.py <- system.file('py', 'utils.py', package = 'chromWaveR')
  source_python(utils.py)
  # One-hot encode DNA strings
  oneHot.strs <- one_hot_encode(as.character(dna.strs))
  return(list(oneHot.strs))
}

#' Function to transform chromWave predictions to nucleosome occupancy.
#'
#' Chromwave models were trained to perform multi class classification,
#' meaning they give n class probabilites back.
#' The function takes the class probabilites as input and transforms them to a
#' continuous signal which can be interpreted as nucleosome occupancy
#' (Nuc/TF occupancy in case of the TF-NUC model).
#' Transformation is performed with the same parameters as used in the chromWave paper.
#'
#' @param \code{preds} Predictions made with a chromWave model.
#' @param \code{model.name} Model name of the chromWave model.
#'
#' @importFrom reticulate source_python
#'
#' @return Function return predicted nucleosome occupancy.
#'
#' @export
#'
#' @examples
#' # Get some DNA sequences and one-hot encode them
#' dna.strs <- rmdDNA()
#' dna.strs <- oneHotEncode(dna.strs)
#'
#' # Load model
#' model <- loadChromWaveModel('invitro')
#'
#' # Perform prediciton
#' dna.preds <- model$predict(dna.strs)
#'
#' # Transform predictions to nuc occ
#' dna.preds <- transformPredictionsToOccupancies(dna.preds, model.name = 'invitro')
#'
transformPredictionsToOccupancies <- function(preds, model.name) {
  # Check for valid model.name
  if(!(model.name %in% availableModels()$model_name)) {
    stop('Invalid model.name.
         Please check availableModels() for valid model names.')
  }

  # Load functions from *.py
  smooth.py <- system.file('py', 'smooth_functions.py', package = 'chromWaveR')
  source_python(smooth.py)

  # Get pre-processing parameters for chromWave model used to derive predictions
  preprocess.params <- get_preprocessing_parameters(model_type = model.name)

  # Transform predicted nuc occupancy
  if(!is.list(preds)) preds <- list(preds)
  nuc.occs <- invert_discretizing(preds, preprocess.params)

  return(nuc.occs)
}
