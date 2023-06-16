#' Plot the Guide Strand with different optional seeds
#'
#' @param guide.seq Guide a.k.a anti-sense sequence oriented 5' > 3'. Sequence must
#' be greater than 8 bp.
#'
#' @return A msaggplot of the guide sequence in addition to the available
#' seed sequences
#'
#' @export
#'
#' @examples
#' library(msa)
#'
#' # Ttr siRNA sequence
#' guide.seq = "UUAUAGAGCAAGAACACUGUUUU"
#'
#' # generate seed plot
#' plotted.seeds = plot_seeds(guide.seq)
plot_seeds <- function(guide.seq){
    mer7A1 = get_seed(guide.seq, "mer7A1")
    mer8 = get_seed(guide.seq, "mer8")
    mer6 = get_seed(guide.seq, "mer6")
    mer7m8 =  get_seed(guide.seq, "mer7m8")

    # Define the set of seed sequences
    set.seq = list("Guide" = mer7A1$Guide,
                 "mer7A1" = mer7A1$Seed.Seq.RNA,
                 "mer8" = mer8$Seed.Seq.RNA,
                 "mer6" = mer6$Seed.Seq.RNA,
                 "mer7m8" = mer7m8$Seed.Seq.RNA)

    # Set the names of the sequence set
    names(set.seq) = c("Guide", "mer7A1", "mer8", "mer6", "mer7m8")

    # Create a RNAStringSet
    plot.seqs.set = Biostrings::RNAStringSet(set.seq, use.names = TRUE)

    # Perform multiple sequence alignment
    results.msa = msa::msa(plot.seqs.set)

    # Get the rna string set from the resultsing search
    rna.str = Biostrings::RNAStringSet(results.msa)

    msa.plot = ggmsa::ggmsa(rna.str,
                          font = 'DroidSansMono',
                          color = "Chemistry_NT",
                          none_bg = FALSE,
                          seq_name = TRUE,
                          use_dot = FALSE,
                          ref = "Guide",
                          consensus_views = TRUE,
                          by_conservation= FALSE) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face="bold",
                                                       color="black",
                                                       size=16),
                   axis.text.x = ggplot2::element_text(size=14)) +
    ggplot2::scale_x_continuous(breaks=seq_len(nchar(guide.seq)))

    return (msa.plot)
}
