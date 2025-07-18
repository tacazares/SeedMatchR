#' Plot a histogram of kmer frequency across genes
#'
#' @param seqs Input [Biostrings::DNAStringSet]
#' @param target.seq Character vector of target seed sequence. Must be same size as width parameter.
#' @param width Width of the kmer to count.
#' @param step Step size for kmer window along the sequence.
#' @param title Figure title
#'
#' @return A list of:
#' *`$plot`: A histogram of frequencies X number of genes with kmer
#' *`$matrix`: Counts matrix
#' @export
#'
#' @examplesIf interactive()
#' annodb = load_annotations("hg38", canonical = T)
#'
#' kmer_hist(annodb$seqs, "AAAAAAA")
kmer_hist <- function(seqs, target.seq, width = 7, step = 1, title = "K-mer gene frequencies"){
  olig.freq.mat =  Biostrings::oligonucleotideFrequency(x = seqs,
                                                        width = width,
                                                        step = step,
                                                        simplify.as = "matrix")

  olig.freq.bin = olig.freq.mat

  olig.freq.bin[olig.freq.bin > 0] <- 1

  gene.counts = colSums(olig.freq.bin)

  target.counts = gene.counts[target.seq][[1]]

  olig.hist = plot(graphics::hist(gene.counts, breaks = 50, xlab = "Number of genes with k-mer", col = "grey90", main = title))

  olig.hist =  olig.hist + graphics::abline(v=target.counts, col="firebrick", lw = 3)

  return(olig.hist)
}


#' Count the k-mer frequency in a set of sequences
#'
#' @param seqs Input [Biostrings::DNAStringSet]
#' @param tx2gene Data frame mapping tx > gene with columsn tx_id + gene_id
#' @param out.path Output directory path
#' @param path.basename Output basename
#' @param width K-mer size
#' @param step step size
#' @param simplify.as output type as either matrix or list
#'
#' @return oligonucleotide frequencies
#' @export
#'
#' @examplesIf interactive()
sequence_kmer_counts <-
  function(seqs,
           tx2gene,
           out.path = NULL,
           path.basename = NULL,
           width = 7,
           step = 1,
           simplify.as = "matrix") {
    # Calculate the oligonucleotide frequency
    olig.freq = Biostrings::oligonucleotideFrequency(seqs,
                                                     width = width,
                                                     step = step,
                                                     simplify.as = simplify.as)

    if (simplify.as == "matrix") {
      olig.freq = as.data.frame(olig.freq)

      rownames(olig.freq) = names(seqs)

      olig.freq = tibble::rownames_to_column(olig.freq, var = "tx_id")

      olig.freq$gene_id = tx2gene$gene_id[match(olig.freq$tx_id, tx2gene$tx_id)]

      olig.freq = olig.freq %>% dplyr::relocate(gene_id, .after = tx_id)
    }

    if (!is.null(out.path)) {
      # Save data as R object and tsv
      save(olig.freq,
           file = file.path(
             out.path,
             paste0(path.basename, "_", width, "mer_counts", ".Rdata")
           ))

      utils::write.table(olig.freq,
                         file = file.path(
                           out.path,
                           paste0(path.basename, "_", width, "mer_counts", ".tsv")
                         ),
                         sep = "\t")
    }

    return(olig.freq)
  }
