#' Plot matches as alignments against target sequence feature
#'
#' @param match.info A row from the `data.frame` of matches produced by `SeedMatchR`
#' @param guide.seq The full siRNA guide sequence
#' @param seqs The feature sequences used in scanning
#' @param aln.name The column name for the `data.frame` reporting alignments
#' @param left.pad The number of nucleotides to the left to show on the feature sequence beyond the length of the guide sequence.
#' @param right.pad The number of nucleotides to the right to show on the feature sequence beyond the length of the guide sequence.
#' @param seq.start.pos The seed definition start position
#' @param seq.stop.pos The seed definition stop position
#' @param output.type Either `data.frame` or `ggplot`
#'
#' @return A `data.frame` representing the input sequences and alignment between them.
#' @export
plotAlignment <- function(match.info,
                          guide.seq,
                          seqs,
                          aln.name = "alignment",
                          left.pad = 0,
                          right.pad = 0,
                          seq.start.pos = 2,
                          seq.stop.pos = 8,
                          output.type = c("data.frame", "ggplot")){
  x <- y <- label <- NULL

  output.type = match.arg(output.type)

  feature.seq = seqs[names(seqs) %in% match.info$names]

  guide.len = nchar(guide.seq)

  full.guide = get_seed(guide.seq, "Full", start.pos = 1, stop.pos = guide.len)

  sirna.target.str = paste0(.make_space(left.pad),
                            as.character(full.guide$Target.Seq),
                            .make_space(right.pad))

  sirna.str.base = as.character(Biostrings::complement(full.guide$Target.Seq))

  d1 = (guide.len - seq.stop.pos)
  d2 = (seq.start.pos - 1)

  matching.seq = XVector::subseq(feature.seq,
                        match.info[1,"start"] - d1 - left.pad,
                        match.info[1,"end"] + d2 + right.pad)

  matching.str = as.character(matching.seq)

  match_char = .matchStrings(sirna.target.str,
                             matching.str)

  aln.df = .create_alignment_df(sirna.str.base,
                                match_char,
                                matching.str,
                                aln.name,
                                left.pad,
                                right.pad)

  if (output.type == "data.frame"){
    return(aln.df)

  } else {
    d2 <- paste0(paste(c("miRNA ","      ","target"), aln.df$alignment), collapse="\n")

    aln.plot <- ggplot2::ggplot(data.frame("x"=1, "y"=1, "label"=d2), ggplot2::aes(x,y,label=label)) +
      ggplot2::theme_void() + ggplot2::geom_text(family="mono", fontface="bold")

    return(aln.plot)

  }
}

#' Match two strings of the same size and return a character vector to represent positions that match
#'
#' @param s1 String 1 to match
#' @param s2 String 2 to match
#'
#' @return A character vector with `|` in the index for every match between both strings.
#' @export
#'
#' @examples
#' SeedMatchR:::.matchStrings("AAAACAGTGTTCTTGCTCTATAA", "ATAACAGTCTTCTTGCTCTGTAA")
.matchStrings <- function(s1, s2){
  stopifnot(nchar(s1) == nchar(s2))

  s1 <- strsplit(s1,"")[[1]]
  s2 <- strsplit(s2,"")[[1]]
  mm <- ifelse(s1==s2, "|", " ")
  mm <- gsub("\\-$", " ", paste(mm, collapse=""))
  mm <- gsub("\\-[^|]|[^|]\\-", "  ", mm)
  mm
}

#' Create a vector of spaces for padding sequences for alignment
#'
#' @param num.spaces Number of spaces to create in vector.
#'
#' @return A character vector of spaces the size of `num.spaces`
#' @export
.make_space <- function(num.spaces){
  paste0(rep(" ",num.spaces),
         collapse="")
}

#' Create a dataframe of alignment information
#'
#' @param sirna The siRNA string as a character vector
#' @param match.string The match characters as a character vector
#' @param feature.seq The feature sequence as a character vector
#' @param aln.name The column name to use for the DF
#' @param left.pad The amount of left spaces to add based on upper-level arguments
#' @param right.pad The amount of right spaces to add based on upper-level arguments
#'
#' @return A `data.frame` that represents the alignment of the siRNA and feature sequence
#' @export
.create_alignment_df <- function(sirna, match.string, feature.seq, aln.name, left.pad, right.pad) {
  stopifnot(is.character(sirna),
              is.character(match.string),
              is.character(feature.seq))

  sirna.str = paste0(.make_space(5 + left.pad),
                       "3'-",
                       sirna,
                       "-5'",
                       .make_space(3 + right.pad))

  match.str = paste0(.make_space(8), match.string, .make_space(6))

  feature.str = paste0("5'-...", feature.seq, "...-3'")

  feature.str = gsub("T", "U", feature.str)
  sirna.str = gsub("T", "U", sirna.str)
  match.str = gsub("T", "U", match.str)

  aln.df <- data.frame(row.names=c("siRNA ",
                                .make_space(6),
                                "target"),
                    "alignment"=c(sirna.str,
                                  match.str,
                                  feature.str))

  aln.df$alignment <- paste0(vapply(max(nchar(aln.df$alignment))-nchar(aln.df$alignment),
                                     FUN.VALUE=character(1),
                                     FUN=function(x) paste0(rep(" ",x),collapse="")),
                             aln.df$alignment)

  return(aln.df)

}
