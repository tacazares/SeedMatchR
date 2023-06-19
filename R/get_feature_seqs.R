#' Get transcripts features and feature sequences
#'
#' @description This function is used to get the genomic features of interest
#' and the DNA sequences associated with them. This function takes advantage of
#' the `GenomicFeatures` package functions `threeUTRsByTranscript`,
#' `fiveUTRsByTranscript`, `exonsBy`, `intronsByTranscript`, and `cdsBy`. These
#' functions are used to generate the features given an input `tx.db` object. A
#' 2bit `dna` input is also required for extracting features sequences.
#'
#' The output of the this function is:
#'
#' * `$db`: the feature `GRanges` object
#' * `$seqs`: `DNAStringSet` of sequences associated to those features
#'
#' @param tx.db A tx.db object
#' @param dna A 2bit dna sequence
#' @param feature.type The type of feature to return. Options: 3UTR, 5UTR,
#' exons, introns, cds
#'
#' @return list containing the feature db object and the feature sequences
#' @export
#' @examplesIf interactive()
#' anno.db = load_species_anno_db("human")
#'
#' features = get_feature_seqs(anno.db$tx.db, anno.db$dna)
get_feature_seqs <- function(tx.db, dna, feature.type="3UTR"){

  # Extract feature based on feature.type argument
  if (feature.type == "3UTR") {
    feature.db = GenomicFeatures::threeUTRsByTranscript(tx.db,
                                                        use.names=TRUE)

  } else if (feature.type == "5UTR") {
    feature.db = GenomicFeatures::fiveUTRsByTranscript(tx.db,use.names=TRUE)

  } else if (feature.type == "exons") {
    feature.db = GenomicFeatures::exonsBy(tx.db, by="tx", use.names=TRUE)

  } else if (feature.type == "introns") {
    feature.db = GenomicFeatures::intronsByTranscript(tx.db, use.names=TRUE)

  } else if (feature.type == "cds"){
    feature.db = GenomicFeatures::cdsBy(tx.db, by="tx", use.names=TRUE)

  } else {
    warning("You have failed to choose a valid option.")
  }

  # Extract the 3p utr sequences for each transcript
  feature.seq = GenomicFeatures::extractTranscriptSeqs(dna, feature.db)

  # Return a list containing the feature database and feature sequences
  results.list = list(db = feature.db,
                      seqs = feature.seq)

  return(results.list)
}


