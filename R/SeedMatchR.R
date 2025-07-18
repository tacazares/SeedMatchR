#' Find seed matches in genomic features
#'
#' @description Find seed matches in a [Biostrings::DNAStringSet] object of
#' sequences.
#'
#' This function will use `get.seed` to extract the seed sequence
#' from the guide sequence.
#'
#' There are two modes for running `SeedMatchR`. The different modes ares set by
#' the `res.format` argument with the options `DESEQ2` or `data.frame`.
#'
#' The `DESEQ2` mode will match the input seed across all rows of the
#' [Biostrings::DNAStringSet] object using [Biostrings::vcountPattern()]. These matches are then
#' aggregated and matched to genes in the `DESEQ2` results `data.frame`. Counts
#' are reported as an additional column with the seed count. By default, this
#' mode will be run if a `DESEQ2` results `data.frame` is provided to the `res` argument.
#'
#' The `data.frame` mode will return a `data.frame` with the interval ranges
#' associated with matches from [Biostrings::vmatchPattern()]. This is the default mode run
#' when no `DESEQ2` results are provided.
#'
#' @inheritParams .get_counts_column
#' @inheritParams .get_match_df
#' @inheritParams get_seed
#' @param res An optional `DESEQ2` results `data.frame`. If provided, an
#' additional column with the seed match count will be added to the
#' `data.frame`. If not provided, `SeedMatchR` will return the interval ranges
#' for each match for the input [Biostrings::DNAStringSet].
#' @param seqs The [Biostrings::DNAStringSet] object with sequence information
#' for features. The names of the sequences should be the transcript names.
#' @param sequence The [Biostrings::DNAString] guide sequence oriented 5' > 3'.
#' @param res.format Format for the returned results. Either 'data.frame', 'DESEQ2', 'iranges', or 'granges'
#' @param start.pos The seed start position
#' @param stop.pos The seed stop position
#' @param match.df Optional: If a matches df is provided the results of the
#' current search will be added with rbind.
#' @param fixed Require that each sequence symbol matches when searching. Should be FALSE if using wobbles.
#' @param get_seed If True, parse the input character vector and return get_seed object.
#' @param max.mismatch Number of allowed mismatches or the total edit distance
#' @param with.indels If True, include indels
#'
#' @return Either a 'data.frame', 'DESEQ2', 'iranges', or 'granges' object based on res.format input
#'
#' @export
#'
#' @examplesIf interactive()
#' library(dplyr)
#'
#' seq = "UUAUAGAGCAAGAACACUGUUUU"
#'
#' anno.db = load_annotations("rnor7")
#'
#' get_example_data("sirna")
#'
#' sirna.data = load_example_data("sirna")
#'
#' res <- sirna.data$Schlegel_2022_Ttr_D1_30mkg
#'
#' # Filter DESeq2 results for SeedMatchR
#' res = filter_res(res, fdr.cutoff=1, fc.cutoff=0, rm.na.log2fc = TRUE)
#'
#' res = SeedMatchR(res = res, seqs = anno.db$seqs,
#' sequence = guide.seq, seed.name = "mer7m8")
SeedMatchR <- function(seqs,
                       sequence,
                       res = NULL,
                       res.format = c("DESEQ2", "data.frame", "granges", "iranges"),
                       seed.name = "mer7m8",
                       col.name = NULL,
                       start.pos = NULL,
                       stop.pos = NULL,
                       match.df = NULL,
                       sirna.name = NULL,
                       shared_genes = TRUE,
                       allow_wobbles = FALSE,
                       get_seed = TRUE,
                       max.mismatch = 0,
                       with.indels = TRUE,
                       fixed = TRUE){
  res.format <- match.arg(res.format)

  if(get_seed){
      # Get the seed sequence of interest
  seed.seq = get_seed(guide.seq = sequence,
                      seed.name = seed.name,
                      start.pos = start.pos,
                      stop.pos = stop.pos,
                      allow_wobbles = allow_wobbles)
  }else{
    seed.seq = sequence
  }


  if (res.format == "DESEQ2" & !(is.null(res))){
    seed.match.res = .get_counts_column(res = res,
                                        target.seq = seed.seq$Target.Seq,
                                        seqs = seqs,
                                        seed.name = seed.name,
                                        col.name = col.name,
                                        shared_genes = shared_genes,
                                        max.mismatch = max.mismatch,
                                        with.indels = with.indels,
                                        fixed = fixed
    )

  } else if (res.format == "DESEQ2" & is.null(res)){
    stop("Result format request is DESEQ2, but no input provided!")

  } else if (res.format %in% c("data.frame", "granges")){
    seed.match.res = .get_match_df(target.seq = seed.seq$Target.Seq,
                                   seqs = seqs,
                                   seed.name = seed.name,
                                   sirna.name = sirna.name,
                                   max.mismatch = max.mismatch,
                                   with.indels = with.indels,
                                   fixed = fixed
    )

    # If a match.df was provided, concat with current results
    if (!is.null(match.df)){
      seed.match.res = rbind(seed.match.res, match.df)
    }

    if (res.format == "granges"){
      seed.match.res = GenomicRanges::GRanges(seed.match.res)
    }

  } else if (res.format == "iranges"){
    seed.match.res = unlist(Biostrings::vmatchPattern(pattern = seed.seq$Target.Seq,
                              subject = seqs,
                              max.mismatch = max.mismatch,
                              with.indels = with.indels,
                              fixed = fixed))

  } else {
    stop("Error!")
  }

  return(seed.match.res)
}


#' Get match results as a column of counts appended to a differential expression results `data.frame`
#'
#' @param res Input differential expression results `data.frame`
#' @param target.seq Pattern to search for in `seqs`
#' @param seqs The `DNAStringSet` to be search for instances of the `target.seq`
#' @param seed.name The seed name to be reported in the `data.frame` column
#' called `seed`
#' @param col.name The optional name of the column of match counts. Will default
#' to `seed.name` if not set
#' @param shared_genes If true, the tx set is reduced to overlapping
#' features in both the sequences db and in the DE data.
#' @param max.mismatch Number of allowed mismatches or the total edit distance
#' @param with.indels If True, include indels
#' @param fixed Require that each sequence symbol matches when searching. Should be FALSE if using wobbles.
#'
#' @return A column of match counts for each gene in `res`
#' @export
#'
#' @examplesIf interactive()
#'
#' guide.seq = "UUAUAGAGCAAGAACACUGUUUU"
#'
#' get_example_data("sirna")
#'
#' sirna.data = load_example_data("sirna")
#'
#' res <- sirna.data$Schlegel_2022_Ttr_D1_30mkg
#'
#' #Load test data
#' res <- Schlegel_2022_Ttr_D1_30mkg
#'
#' # Filter results for SeedMatchR
#' res = filter_res(res, fdr.cutoff=1, fc.cutoff=0, rm.na.log2fc = TRUE)
#'
#' seed.seq = get_seed(guide.seq = guide.seq, seed.name = "mer7m8")
#'
#' anno.db = load_annotations("rnor7")
#'
#' SeedMatchR:::.get_counts_column(res = res, target.seq = seed.seq$Target.Seq,
#' seqs = anno.db$seqs, seed.name = "mer7m8")
.get_counts_column <- function(res,
                               target.seq,
                               seqs,
                               seed.name,
                               col.name = NULL,
                               shared_genes = TRUE,
                               max.mismatch = 0,
                               with.indels = F,
                               fixed = F){
  if (shared_genes == TRUE){
    res = res %>% dplyr::filter(gene_id %in% names(seqs))
  }

  seed.matches = Biostrings::vcountPattern(pattern = target.seq,
                                           subject = seqs,
                                           max.mismatch = max.mismatch,
                                           with.indels = with.indels,
                                           fixed = fixed)

  # Create a data.frame of match counts per transcript
  match.df = data.frame(gene_id = names(seqs),
                        match.counts = seed.matches)

  # Match the rows in the summary df to deseq.genes and report the counts back
  match.counts <- match.df$match.counts[match(res$gene_id,
                                              match.df$gene_id)]

  # fill any missing na values with 0
  match.counts[is.na(match.counts)] = 0

  # Concat the seed match results
  res = cbind(res, seed.name = match.counts)

  names(res)[names(res) == 'seed.name'] <- if (is.null(col.name))
                                               seed.name else col.name

  return (res)
}


#' Get sequence matches as a `data.frame`
#'
#' This function will perform the siRNA sequence match to the transcript library
#' and return it as a `data.frame` of ranges. This function uses
#' `vmatchPattern` to perform the matching. The `data.frame` will have
#' the columns: "names", "start", "end", "width", "seed".  If the optional
#' `sirna.name` argument is provided then an additional column will be added.
#'
#' @param target.seq Pattern to search for in `seqs`
#' @param seqs The `DNAStringSet` to be search for instances of the `target.seq`
#' @param seed.name The seed name to be reported in the `data.frame` column called `seed`
#' @param sirna.name Optional siRNA name. A new column called `sirna.name` will be added to `data.frame`
#' @param max.mismatch Number of allowed mismatches or the total edit distance
#' @param with.indels If True, include indels
#' @param fixed Require that each sequence symbol matches when searching. Should be FALSE if using wobbles.
#'
#' @return A `data.frame` of matches and ranges for each sequence in `seqs`
#' @export
#'
#' @examplesIf interactive()
#'
#' guide.seq = "UUAUAGAGCAAGAACACUGUUUU"
#'
#' seed.seq = get_seed(guide.seq = guide.seq, seed.name = "mer7m8")
#'
#' anno.db = load_annotations("rnor7")
#'
#' SeedMatchR:::.get_match_df(target.seq = seed.seq$Target.Seq,
#' seqs = anno.db$seqs, seed.name = "mer7m8")
.get_match_df <- function(target.seq,
                          seqs,
                          seed.name,
                          sirna.name = NULL,
                          max.mismatch = 0,
                          with.indels = TRUE,
                          fixed = TRUE){

  seq.matches = vmatchPattern2(
    pattern = target.seq,
    subject = seqs,
    max.mismatch = max.mismatch,
    with.indels = with.indels,
    fixed = fixed
  )

  seq.matches = seq.matches %>% dplyr::rename(seqnames = names)

  seq.matches = seq.matches[c("seqnames", "start", "end", "width")]

  if (!is.null(sirna.name)){
    seq.matches$sirna = sirna.name
    seq.matches$query = seed.name
    seq.matches$TargetSeq = as.character(target.seq)
    #seq.matches$MatchString = apply(seq.matches, 1, get_match_seq, seqs = seqs)
    seq.matches$EditDistance = max.mismatch
    seq.matches$AllowIndels = with.indels
  }

  return(seq.matches)

}
