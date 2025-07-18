utils::globalVariables(c("mm","query"))

#' Search input `Biostrings::DNAStringSet` for a set of defined sequences
#'
#' This function will use `Biostrings::vcountpattern` to find matches to the
#' input siRNA sequence. This function will then provide a matrix that is
#' sequences x match type. The best match will be reported based on a pre-defined
#' list of rankings.
#'
#' @param guide.seq Input guide sequence.
#' @param seqs `Biostrings::DNAStringSet` of sequences to search.
#' @param group.name Group name to use for labeling.
#' @param indel.bool Boolean in long format (TRUE, FALSE) for whether to include indels.
#' @param seqs.to.search List of sequences to search. Options: "Full", "18mer", "15mer", "mer8", "mer7m8", "mer7A1", "mer6".
#' @param res_format Either "iranges" or "counts".
#' @param allow_wobbles If true, all positions for a G:U wobble are considered.
#' @param fixed Whether to use matching string symbols for search when ambiguous.
#' @param nonseed.edit.distances An integer value indicating allowed edit distance for non-seed searches.
#' @param seed.edit.distances An integer value indicating allowed edit distance for seed searches.
#'
#' @return A `data.frame` of sequence match counts.
#' @export
#'
#' @examplesIf interactive()
#' print("TBD")
full_search <- function(guide.seq,
                        seqs,
                        group.name = "siRNA",
                        indel.bool = FALSE,
                        allow_wobbles = FALSE,
                        seqs.to.search = c("Full", "18mer", "15mer", "mer8", "mer7m8", "mer7A1", "mer6"),
                        res_format = c("counts", "iranges"),
                        fixed = FALSE,
                        nonseed.edit.distances = c(0:4),
                        seed.edit.distances = c(0:1)){
  res_format <- match.arg(res_format)
  stopifnot(sum(seqs.to.search %in% c("Full", "18mer",  "15mer", "mer8", "mer7m8", "mer7A1", "mer6")) > 0)

  match_matrix = data.frame(seqs = names(seqs))
  match.df = c()

  for (si.seq in seqs.to.search){

    if (si.seq %in% c("Full", "18mer", "15mer")){
      edit.dist.list = nonseed.edit.distances
      with.indel = indel.bool
    } else {
      edit.dist.list = seed.edit.distances
      with.indel = indel.bool
    }

    as.seq = get_seed(guide.seq, seed.name = si.seq, allow_wobbles = allow_wobbles)

    for (edit.dist in edit.dist.list){
      search.name = paste0(si.seq, ".mm", edit.dist)

      if (res_format == "counts"){
        matches = Biostrings::vcountPattern(as.seq$Target.Seq,
                                               seqs,
                                               max.mismatch = edit.dist,
                                               with.indel = with.indel,
                                               fixed = fixed)

        match_matrix = cbind(match_matrix, matches)

        names(match_matrix)[ncol(match_matrix)] <- search.name

        match.df <- match_matrix

      } else {
        matches = vmatchPattern2(pattern = as.seq$Target.Seq,
                                  subject = seqs,
                                  max.mismatch = edit.dist,
                                  with.indels = with.indel,
                                  fixed = fixed)

        if (nrow(matches) > 0) {
          matches$query = search.name
        }

        match.df = rbind(match.df, matches)

      }
    }
  }

  match.df$Group = group.name

  if (res_format == "counts"){
    match.df = .get_best_match(match.df)
  }

  return(match.df)
}

.get_best_match <- function(match_matrix){
  #TODO find a better way.
  match_matrix = match_matrix %>%
      dplyr::mutate(
        best.match = dplyr::case_when(
          Full.mm0 >= 1 ~ "Full-D0",
          Full.mm1 >= 1 ~ "Full-D1",
          Full.mm2 >= 1 ~ "Full-D2",
          Full.mm3 >= 1 ~ "Full-D3",
          Full.mm4 >= 1 ~ "Full-D4",
          `18mer.mm0` >= 1 ~ "18mer-D0",
          `18mer.mm1` >= 1 ~ "18mer-D1",
          `18mer.mm2` >= 1 ~ "18mer-D2",
          `18mer.mm3` >= 1 ~ "18mer-D3",
          `18mer.mm4` >= 1 ~ "18mer-D4",
          `15mer.mm0` >= 1 ~ "15mer-D0",
          `15mer.mm1` >= 1 ~ "15mer-D1",
          `15mer.mm2` >= 1 ~ "15mer-D2",
          `15mer.mm3` >= 1 ~ "15mer-D3",
          `15mer.mm4` >= 1 ~ "15mer-D4",
          mer8.mm0 >= 1 ~ "8mer",
          mer7m8.mm0 >= 1 ~ "7mer-m8",
          mer7A1.mm0 >= 1 ~ "7mer-A1",
          mer6.mm0 >= 1 ~ "6mer",
          .default = "No.matches"
        )
      )

  return(match_matrix)
}

#' Search input `Biostrings::DNAStringSet` for a set of defined sequences
#'
#' This function will use `Biostrings::vcountpattern` to find matches to the
#' input siRNA sequence. This function will then provide a matrix that is
#' sequences x match type. The best match will be reported based on a pre-defined
#' list of rankings.
#'
#' @param guide.seq Input guide sequence
#' @param seqs `Biostrings::DNAStringSet` of sequences to search
#' @param group.name Group name to use for labeling
#' @param indel.bool.opt Vector of indel options.
#' @param seqs.to.search List of sequences to search. Options: "Full", "18mer", "mer8", "mer7m8", "mer7A1", "mer6"
#' @param edit.dist.list A vector of possible edit distances.
#' @param fixed.opt A vector of possible fixed paramters.
#' @param max.seed.edit.dist An integer value indicating allowed edit distance for seed searches.
#'
#' @return A `data.frame` of sequence match counts
#' @export
#'
#' @examplesIf interactive()
#' print("TBD")
full_search_ranges <- function(guide.seq,
                               seqs,
                               group.name = "siRNA",
                               indel.bool.opt = c(TRUE, FALSE),
                               edit.dist.list = c(0:4),
                               seqs.to.search = c("Full", "18mer", "15mer", "mer8", "mer7m8", "mer7A1", "mer6"),
                               fixed.opt = FALSE,
                               max.seed.edit.dist = 1) {
  df <- expand.grid(indel.bool.opt,
                    seqs.to.search,
                    edit.dist.list)

  colnames(df) <- c("indels", "query", "mm")

  df = df %>% dplyr::filter(!(query %in% c("mer8", "mer7m8", "mer7A1", "mer6") & mm > max.seed.edit.dist))

  df$query = as.character(df$query)

  seq.library = stats::setNames(lapply(seqs.to.search, get_seed, guide.seq = guide.seq), seqs.to.search)

  seq.library = sapply(seq.library, function(x) x[["Target.Seq"]])

  out.res = mapply(
    .run_scanning,
    si.seq = df$query,
    edit.dist = df$mm,
    with.indel = df$indels,
    MoreArgs = list(
      seq.library = seq.library,
      seqs = seqs,
      fixed = fixed.opt
    ))

  out.res = lapply(out.res, function(df) stats::na.omit(df))

  out.df = dplyr::bind_rows(out.res)

  out.df$Group = group.name

  return(out.df)
}


.run_scanning <- function(seq.library,
                          si.seq,
                          seqs,
                          edit.dist,
                          with.indel,
                          fixed) {
  as.seq = seq.library[[si.seq]]

  search.name = paste0(si.seq,
                       ".mm",
                       edit.dist,
                       ".indel",
                       with.indel)

  out.matches = vmatchPattern2(
    pattern = as.seq,
    subject = seqs,
    max.mismatch = edit.dist,
    with.indels = with.indel,
    fixed = fixed
  )

  if (nrow(out.matches) > 0) {
    out.matches$query = search.name
    out.matches$TargetSeq = as.character(as.seq)
    out.matches$seqnames = out.matches$names
    #out.matches$MatchString = apply(out.matches, 1, get_match_seq, seqs = seqs)
    out.matches$EditDistance = edit.dist
    out.matches$AllowIndels = with.indel
  }

  return(out.matches)
}

.get_match_seq = function(row, seqs){
  tx.id = row["seqnames"]
  str.match = substr(seqs[tx.id], row["start"], row["end"])
  return(str.match)
}

.get_alignment = function(row){
  s1 = row["TargetSeq"]
  s2 = row["MatchString"]

  alg <- Biostrings::pairwiseAlignment(s1, s2)

  return(as.character(Biostrings::alignedSubject(alg)))
}

