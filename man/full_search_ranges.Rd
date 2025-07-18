% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/full_search.R
\name{full_search_ranges}
\alias{full_search_ranges}
\title{Search input \code{Biostrings::DNAStringSet} for a set of defined sequences}
\usage{
full_search_ranges(
  guide.seq,
  seqs,
  group.name = "siRNA",
  indel.bool.opt = c(TRUE, FALSE),
  edit.dist.list = c(0:4),
  seqs.to.search = c("Full", "18mer", "15mer", "mer8", "mer7m8", "mer7A1", "mer6"),
  fixed.opt = FALSE,
  max.seed.edit.dist = 1
)
}
\arguments{
\item{guide.seq}{Input guide sequence}

\item{seqs}{\code{Biostrings::DNAStringSet} of sequences to search}

\item{group.name}{Group name to use for labeling}

\item{indel.bool.opt}{Vector of indel options.}

\item{edit.dist.list}{A vector of possible edit distances.}

\item{seqs.to.search}{List of sequences to search. Options: "Full", "18mer", "mer8", "mer7m8", "mer7A1", "mer6"}

\item{fixed.opt}{A vector of possible fixed paramters.}

\item{max.seed.edit.dist}{An integer value indicating allowed edit distance for seed searches.}
}
\value{
A \code{data.frame} of sequence match counts
}
\description{
This function will use \code{Biostrings::vcountpattern} to find matches to the
input siRNA sequence. This function will then provide a matrix that is
sequences x match type. The best match will be reported based on a pre-defined
list of rankings.
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
print("TBD")
\dontshow{\}) # examplesIf}
}
