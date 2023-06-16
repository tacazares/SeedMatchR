#' Plot a histogram of dstat statistics from tests provided by the `twosamples`
#' package
#'
#' @description This function will take the twosamples object and plot the
#' histogram of dstats using ggplot.
#'
#' @param test.results The twosamples object holding the test stats and
#' iteration information.
#' @param title Title of the plot
#' @param line.col Color to use for Dstat line
#'
#' @return ggplot histogram object
#' @export
#' @examplesIf interactive()
#' library(dplyr)
#'
#' guide.seq = "UUAUAGAGCAAGAACACUGUUUU"
#'
#' anno.db = load_species_anno_db("human")
#'
#' features = get_feature_seqs(anno.db$tx.db, anno.db$dna)
#'
#' # Load test data
#' get_example_data("sirna")
#'
#' load_example_data("sirna")
#'
#' # Filter DESeq2 results for SeedMatchR
#' res = filter_deseq(res, fdr.cutoff=1, fc.cutoff=0, rm.na.log2fc = TRUE)
#'
#' res = SeedMatchR(res, anno.db$gtf, features$seqs, guide.seq, "mer7m8")
#'
#' # Gene set 1
#' mer7m8.list = res$gene_id[res$mer7m8 >= 1]
#'
#' # Gene set 2
#' background.list = res$gene_id[!(res$mer7m8 %in% mer7m8.list)]
#'
#' ecdf.res = ecdf_stat_test(res, mer7m8.list, background.list)
#'
#' hist = dstat_histogram(ecdf.res)
dstat_histogram <- function(test.results,
                            title = "Histogram of Test Statistics",
                            line.col = "red"){
  boostraps <- NULL
  # Extract the bootstrap values from the twosamples class
  boostrap_vals = unlist(attributes(test.results)["bootstraps"])

  # Create a histogram df
  hist.df = data.frame("boostraps" = boostrap_vals)

  # Plot annotated histograms
  ggplot.hist = ggplot2::ggplot(hist.df, ggplot2::aes(x=boostraps)) +
    ggplot2::geom_histogram(alpha=.8, colour="black", fill = "grey") +
    ggplot2::geom_vline(ggplot2::aes(xintercept=test.results["Test Stat"]),
                        color=line.col, linetype="dashed", size=1) +
    ggplot2::labs(title = title,
                  x="Dstat",
                  y="Count",
                  subtitle = paste0("Test: ",
                                    attributes(test.results)["test_type"],
                                    "\n",
                                    "Bootstraps: ",
                                    length(boostrap_vals),
                                    "\t",
                                    "Dstat: ",
                                    round(test.results["Test Stat"], 4),
                                    "\t",
                                    "P-Value: ",
                                    round(test.results["P-Value"], 4))) +
    SeedMatchR.theme

  return(ggplot.hist)
}
