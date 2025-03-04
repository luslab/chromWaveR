% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{transformPredictionsToOccupancies}
\alias{transformPredictionsToOccupancies}
\title{Function to transform chromWave predictions to nucleosome occupancy.}
\usage{
transformPredictionsToOccupancies(preds, model.name)
}
\arguments{
\item{\code{preds}}{Predictions made with a chromWave model.}

\item{\code{model.name}}{Model name of the chromWave model.}
}
\value{
Function return predicted nucleosome occupancy.
}
\description{
Chromwave models were trained to perform multi class classification,
meaning they give n class probabilites back.
The function takes the class probabilites as input and transforms them to a
continuous signal which can be interpreted as nucleosome occupancy
(Nuc/TF occupancy in case of the TF-NUC model).
Transformation is performed with the same parameters as used in the chromWave paper.
}
\examples{
# Get some DNA sequences and one-hot encode them
dna.strs <- rmdDNA()
dna.strs <- oneHotEncode(dna.strs)

# Load model
model <- loadChromWaveModel('invitro')

# Perform prediciton
dna.preds <- model$predict(dna.strs)

# Transform predictions to nuc occ
dna.preds <- transformPredictionsToOccupancies(dna.preds, model.name = 'invitro')

}
