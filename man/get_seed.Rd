% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_seed.R
\name{get_seed}
\alias{get_seed}
\title{Get the target sequence for a given siRNA input sequence}
\usage{
get_seed(
  guide.seq,
  seed.name = "mer7m8",
  start.pos = NULL,
  stop.pos = NULL,
  allow_wobbles = FALSE
)
}
\arguments{
\item{guide.seq}{A character string greater than 8 bp and oriented 5'-> 3'.}

\item{seed.name}{The seed name of interest. Options: mer8,  mer7A1, mer7m8,
mer6, 18-mer,  15-mer, or Full. If not in the default list, the start.pos and stop.pos arguments
will be used to define the seed.}

\item{start.pos}{The start position for a custom seed definition. 1-based}

\item{stop.pos}{The stop position for a custom seed definition.}

\item{allow_wobbles}{If True, allow G:U wobbles by replacing U with Y.}
}
\value{
An object with the entries:
\itemize{
\item \code{Guide}: Input guide sequence. Input is expected to be RNA.
\item \code{Seed.Name}: The seed name.
\item \code{Seed.Seq.RNA}: The seed sequence as a \code{RNAString}
\item \code{Seed.Seq.DNA}: The seed sequence as a \code{DNAString}
\item \code{Target.Seq}: The target DNA sequence based on the reverse complement of
the \code{Seed.Seq.DNA}
}
}
\description{
Given a sequence greater than 8 bp oriented 5' -> 3' and a seed
definition (or custom definition), this function will return an object
containing target sequence information.

Built-in options: \code{mer8},  \code{mer7A1}, \code{mer7m8}, \code{mer6}, \code{Full}, and \verb{18mer}

Custom Definitions: Users can input a custom seed name, but must provide the
start position (\code{start.pos}) and stop position (\code{stop.pos}) that define the
range of the seed sequence.

Note: \code{mer8} and \code{mer7A1} force a U at position g1. This results in an A in the target
sequence being searched. This is based on historical data and observations
from miRNA biology. If you do not wish to enforce this, provide the start,
stop, and custom name to use the original sequence.
}
\examples{
# Example Ttr from Schlegel et al. 2022
guide.seq = "UUAUAGAGCAAGAACACUGUUUU"

# Get seed match
seed.seq = get_seed(guide.seq = guide.seq, seed.name = "mer7m8")
}
