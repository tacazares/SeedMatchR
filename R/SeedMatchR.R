#' Find seed matches in genomic features
#'
#' @description Find seed matches in a `DNAStringSet` object of
#' sequences. This function will use `get.seed` extract the seed sequence
#' from the guide sequence. The seed is then searched across all rows of the
#' `DNAStringSet` object using `vpatterncount`.
#'
#' This function returns the input DESeq2 results `data.frame` with an
#' additional column that contains the counts for the input `seed.name`.
#'
#' @param res A DESeq2 results `data.frame`
#' @param gtf GTF file used to map features to genes. The object must have
#' columns transcript_id and gene_id
#' @param seqs The `DNAStringSet` object with sequence information
#' for features. The names of the sequences should be the transcript names.
#' @param sequence The `DNAString` guide sequence oriented 5' > 3'.
#' @param seed.name The name of specific seed to extract. Options are:  mer8,
#' mer7A1, mer7m8, mer6
#' @param col.name The string to use for the column name. Defaults to seed name
#' @param tx.id.col Use the transcript_id column instead of gene_id
#' @param mismatches The number of mismatches to allow in search
#' @param indels Whether to allow indels in search
#'
#' @return A modified DESeq2 results dataframe that has column named after the
#' seed of choice representing the number of match counts.
#'
#' @export
#'
#' @examplesIf interactive()
#' library(dplyr)
#'
#' seq = "UUAUAGAGCAAGAACACUGUUUU"
#'
#' anno.db = load_species_anno_db("human")
#'
#' features = get_feature_seqs(anno.db$tx.db, anno.db$dna)
#'
#' # Load test data
#' res <- Schlegel_2022_Ttr_D1_30mkg
#'
#' # Filter DESeq2 results for SeedMatchR
#' res = filter_deseq(res, fdr.cutoff=1, fc.cutoff=0, rm.na.log2fc = TRUE)
#'
#' res = SeedMatchR(res, anno.db$gtf, features$seqs, seq, "mer7m8")
SeedMatchR <- function(res,
                       gtf,
                       seqs,
                       sequence,
                       seed.name="mer7m8",
                       col.name = NULL,
                       mismatches=0,
                       indels=FALSE,
                       tx.id.col = TRUE){
  gene_id <- NULL
  # Get the seed sequence of interest
  seed.seq = get_seed(sequence, seed.name)

  seed.matches = Biostrings::vcountPattern(seed.seq$Target.Seq,
                                           seqs,
                                           max.mismatch=mismatches,
                                           with.indel=indels)

  if (tx.id.col == TRUE){
    gene_id_names = gtf$gene_id[match(names(seqs), gtf$transcript_id)]
  } else {
    gene_id_names = names(seqs)
  }

  # Create a data.frame of match counts per transcript
  match.df = data.frame(gene_id = gene_id_names,
                        transcript_id = names(seqs),
                        match.counts = seed.matches)

  # Sum the counts by gene_id
  sum_gene_counts <- match.df %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(match_counts = sum(match.counts))

  # Match the rows in the summary df to deseq.genes and report the counts back
  match.counts <- sum_gene_counts$match_counts[match(res$gene_id,
                                                     sum_gene_counts$gene_id)]

  # fill any missing na values with 0
  match.counts[is.na(match.counts)] = 0

  # Concat the seed match results
  res = cbind(res, seed.name = match.counts)

  ifelse(!is.null(col.name),
         names(res)[names(res) == 'seed.name'] <- col.name,
         names(res)[names(res) == 'seed.name'] <- seed.name)

  return (res)
}
