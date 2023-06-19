#' Load species specific \code{AnnotationDb}
#'
#' @description Use `AnnotationHub` to load species-specific GTF and 2bit DNA
#' sequences. This function currently works for human, rat, and mouse.
#'
#' The function will return:
#'
#' * `$gtf`: A `GRanges` object containing the GTF information
#' * `$tx.db`: A `tx.db` object made from the GTF
#' * `$dna`: The 2bit DNA sequence as a `DNAStringSet`
#'
#' @param species.name Species name. Options: human, rat, mouse
#' @param remove.na.rows Remove rows with NA in the gene_id column
#'
#' @return Species specific \code{AnnotationDb}
#' @export
#'
#' @examplesIf interactive()
#' anno.db = load_species_anno_db("human")
load_species_anno_db <- function(species.name, remove.na.rows = TRUE){
  # Load annotation hub
  ah <- AnnotationHub::AnnotationHub()

  message(paste0("Loading annotation data for ", species.name))
  
  if (species.name == "human") {
    gtf <- ah[["AH98495"]]

    dna <- ah[["AH99541"]]

  } else if (species.name == "rat") {
    gtf <- ah[["AH92420"]]

    dna <- ah[["AH93578"]]

  } else if (species.name == "mouse") {
    gtf <- ah[["AH7567"]]

    dna <- ah[["AH49773"]]

  } else {
    stop("Currently, species can only be human, rat, or mouse")
  }

  # filter the gtf for only standard chromosomes
  gtf <- GenomeInfoDb::keepStandardChromosomes(gtf, pruning.mode="coarse")

  ## Drop 'MT'.
  gtf <- GenomeInfoDb::dropSeqlevels(gtf, "MT", pruning.mode="coarse")

  if(remove.na.rows == TRUE){
    gtf = gtf[!is.na(gtf$transcript_id)]
  }

  # Make a tx database
  tx.db = GenomicFeatures::makeTxDbFromGRanges(gtf)

  # Create a list of results to return
  results.list = list(gtf = gtf,
                      tx.db = tx.db,
                      dna = dna)

  return(results.list)
}
