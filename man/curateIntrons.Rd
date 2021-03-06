% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Estim_curateIntrons.R
\name{curateIntrons}
\alias{curateIntrons}
\title{Internal function for excluding low mappability regions from intronic intervals}
\usage{
curateIntrons(
  mappabilityFile = NULL,
  chainPath = NULL,
  readLength,
  minMapScore = 0.5,
  saveDir,
  verbose = TRUE
)
}
\arguments{
\item{mappabilityFile}{Path to the (bed?) file containing mappability score to use.
If NULL, no mappability curation will be performed (Default: NULL).}

\item{chainPath}{(Default: NULL). A \code{Chain} object, usually imported with ['import.chain'], or something coercible to one.
Specify only if genome version liftover of mappability scores is needed.}

\item{readLength}{Consider only low mappability region longer than this value.}

\item{minMapScore}{Masq regions with mappability score below this value (Default: 0.5)}

\item{saveDir}{Path to the directory where all sample results are saved.}

\item{verbose}{Whether information messages should be displayed (DEFAULT: TRUE).}
}
\value{
A \code{GRangesList} having same length as \code{indIntrons} listing sub-intervals
        with mappability score and width above the specified thresholds.
}
\description{
This function removes low mappability regions from intronic intervals.
Blacklisted regions are those having mappability score below a given threshold and spanning a minimal width (eg: read length).
In certain cases (depending on read length and your reference genome) pre-computed mappability scores can be downloaded
from ENCODE databases.
}
