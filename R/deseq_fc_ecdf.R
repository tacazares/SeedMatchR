#' Plot the ecdf for DESeq2 log2 Fold Changes
#'
#' @description This functions will take DESeq2 results as a `data.frame` and
#' plot the ecdf for the input `gene.lists`.
#'
#' The gene sets to plot should be provided as a list of lists.
#'
#' Example:
#'
#' `gene.lists = list("Background" = c("gene1", "gene2"), "Target" = c("gene2",
#' "gene3"), "Overlap" = c("gene2"))`
#'
#' This function will also perform statistical testing using the the
#' `twosamples` package if `plot.hist` is TRUE. The output will be saved to a
#' PDF if an `output.filename` is provided.
#'
#' Users can define the groups that are to be compared in the statistical test
#' using the `null.name` and `target.name` arguments. The names must be found
#' in `gene.lists`. The `factor.order` is used to order the groups in the
#' analysis.
#'
#' This functions returns:
#'
#' * `$plot`: The ECDF plot
#' * `$stats`: The twosamples results object
#'
#' @param res The DESeq2 results dataframe
#' @param gene.lists A nest list of gene names. Example:
#' gene.lists = list("Background" = gene.list2, "Target" = gene.list1,
#' "Overlap" = gene.list3)
#' @param title The tile of the plot
#' @param output.filename If the output filename is provided, then the plot is
#' saved.
#' @param palette The color palette to use for your curves
#' @param factor.order The order to use for the legends
#' @param x.lims The xlimits range
#' @param stats.test The statistic test to use. Options: KS, Kuiper, DTS, CVM,
#' AD, Wass
#' @param nboots The number of iterations to run the simulation
#' @param null.name The name in the gene.list to use as the null for ecdf plots
#' @param target.name The name in the gene.list to use as the target for ecdf
#' @param keep.boots Whether to keep the results from the boot strapping
#' analysis
#' @param plot.hist Whether to plot the histogram of dstats derived from boost
#' strap analysis. Keep.boots must be TRUE.
#' @param height Plot height in inches
#' @param width Plot width in inches
#' @param dpi The dpi resolution for the figure
#'
#' @return A ggplot object of the ecdf
#' @export
#'
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
#' ecdf.results = deseq_fc_ecdf(res,
#' list("Background" = background.list, "mer7m8" = mer7m8.list),
#' stats.test = "KS",
#' factor.order = c("Background", "mer7m8"),
#' null.name = "Background",
#' target.name = "mer7m8",
#' nboots = 100)
deseq_fc_ecdf <- function(res,
                          gene.lists,
                          title="ECDF",
                          output.filename=NULL,
                          palette = SeedMatchR.palette,
                          factor.order=NULL,
                          x.lims = c(-1,1),
                          stats.test = NULL,
                          nboots = 10000,
                          null.name = 1,
                          target.name = 2,
                          keep.boots = TRUE,
                          plot.hist = TRUE,
                          height = 5,
                          width=10,
                          dpi=320){
  # This is to prevent a no visible binding for global variable error
  log2FoldChange <- Group <- NULL

  # Check if gene lists overlap, if true issue warning
  check_gene_list_overlap(gene.lists)

  res$Group  = "NA"

  # For every group name in names(gene.list), label rows based on  list name
  for (group in names(gene.lists)){
    res$Group  = ifelse(res$gene_id %in% unlist(gene.lists[group]),
                        group,
                        res$Group)
  }

  # If the factor order is null, then use the order of the gene.lists names
  # The null.name and target.name should be the first two in the list.
  if(is.null(factor.order)){
    factor.order = names(gene.lists)
  }

  res$Group = factor(res$Group, levels=factor.order)

  # Get the counts of genes in each gene lists
  gene.counts = unlist(lapply(gene.lists, length))

  # Create legend labels with counts
  legend.labels = paste0(factor.order, ": ", gene.counts)

  print(paste0("Comparing: ", null.name, " vs. ", target.name))

  # Use a stats test to determine significance
  test.results = ecdf_stat_test(res,
                                gene.lists[null.name],
                                gene.lists[target.name],
                                stats.test = stats.test,
                                nboots = nboots,
                                keep.boots = keep.boots)

  # Create ggplot
  ecdf.plot = ggplot2::ggplot(as.data.frame(res),
                              ggplot2::aes(x=log2FoldChange,
                                           group=Group,
                                           col=Group,
                                           geom = "step")) +
    ggplot2::stat_ecdf(size = 2) +
    ggplot2::scale_colour_manual(values=palette, labels = legend.labels) +
    ggplot2::labs(title = title,
                  x = "Log2(Fold Change)",
                  y="Cumulative Distribution",
                  subtitle = ifelse(!is.null(stats.test),
                                    paste0("Dstat: ",
                                           round(test.results["Test Stat"],
                                                 3),
                                           "\n",
                                           "Pvalue: ",
                                           test.results["P-Value"]),
                                    "No test stats")) +
    SeedMatchR.theme +
    ggplot2::coord_cartesian(xlim = x.lims)

  # If plot.hist is true, create a two column plot with the ecdf + hist
  if (plot.hist == TRUE){
    testit::assert("If plot.hist == TRUE,then keep.boots must also be TRUE",
                   keep.boots == TRUE)
    hist.plot = dstat_histogram(test.results)
    #final.plot <- ecdf.plot + hist.plot
    final.plot <- cowplot::plot_grid(ecdf.plot, hist.plot)

  } else {
    final.plot <- ecdf.plot
  }

  # If there is an output name provided, save the figure
  if (!is.null(output.filename)){
    ggplot2::ggsave(output.filename,
                    final.plot,
                    device="pdf",
                    height=height,
                    width=width,
                    units = "in",
                    dpi=dpi)
  }

  return(list(plot = final.plot, stats = test.results))
}
