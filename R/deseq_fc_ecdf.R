utils::globalVariables(c("log2FoldChange", "Group"))

#' Plot the ECDF of log2(fold change) between groups
#'
#' @description This functions will take differential expression results as a
#' `data.frame` and plot the ECDF for the input `gene_lists`.
#'
#' The gene sets to plot should be provided as a named list of vectors.
#'
#' Example:
#'
#' `gene_lists = list("Background" = c("gene1", "gene2"), "Target" = c("gene2",
#' "gene3"))`
#'
#' The output will be saved to a PDF if an `output_filename` is provided.
#'
#' Users can define the groups that are to be compared in the statistical test
#' using the `null_name` and `target_name` arguments. The names must be found
#' in `gene_lists`. The `factor_order` is used to order the groups in the
#' analysis.
#'
#' This functions returns:
#'
#' * `$plot`: The ECDF plot
#' * `$stats`: The results object
#'
#' @param res The differential expression results data frame.
#' @param gene_lists A nested list of gene names. Example:
#' gene_lists = list("Background" = background_list, "Target" = target_list)
#' @param title The tile of the plot.
#' @param output_filename If the output filename is provided, then the plot is
#' saved.
#' @param palette Vector of colors to use for your curves.
#' @param factor_order Vector defining the order to use for the legends.
#' @param x_lims The x-axis range limits.
#' @param stats_test The statistical test to use. Options: KS, Wass, Wilcoxen, or
#' DTS.
#' @param alternative The alternative hypothesis to test. Options: greater,
#' less, two.sided.
#' @param null_name The name in the gene.list to use as the null for ECDF plots.
#' @param target_name The name in the gene.list to use as the target for ECDF.
#' @param height Plot height in inches.
#' @param width Plot width in inches.
#' @param dpi The DPI resolution for the figure.
#' @param theme A ggplot2 theme defining plot parameters.
#'
#' @return A ggplot object and ECDF statistics.
#' @export
#'
#' @examplesIf interactive()
#' library(dplyr)
#'
#' guide.seq = "UUAUAGAGCAAGAACACUGUUUU"
#'
#' anno.db = load_annotations("rnor7")
#'
#' # Load test data
#' get_example_data("sirna")
#'
#' sirna.data = load_example_data("sirna")
#'
#' res <- sirna.data$Schlegel_2022_Ttr_D1_30mkg
#'
#' # Filter DESeq2 results for SeedMatchR
#' res = filter_res(res, fdr.cutoff=1, fc.cutoff=0, rm.na.log2fc = TRUE)
#'
#' res = SeedMatchR(res = res, gtf = anno.db$gtf, seqs = anno.db$seqs,
#' sequence = guide.seq, seed.name = "mer7m8", tx.id.col= FALSE)
#'
#' # Gene set 1
#' mer7m8.list = res$gene_id[res$mer7m8 >= 1]
#'
#' # Gene set 2
#' background.list = res$gene_id[res$mer7m8 == 0]
#'
#' ecdf.results = deseq_fc_ecdf(res,
#' list("Background" = background.list, "mer7m8" = mer7m8.list))
deseq_fc_ecdf <- function(res,
                          gene_lists,
                          title = "ECDF",
                          output_filename = NULL,
                          palette = SeedMatchR.palette,
                          factor_order = NULL,
                          x_lims = c(-1, 1),
                          stats_test = c("KS", "Wilcoxen", "Wass", "DTS"),
                          alternative = c("greater", "less", "two.sided"),
                          null_name = 1,
                          target_name = 2,
                          height = 5,
                          width = 5,
                          dpi = 320,
                          theme = SeedMatchR.theme) {
  alternative <- match.arg(alternative)
  stats_test <- match.arg(stats_test)

  # Check if gene lists overlap, if true issue warning
  .check_gene_list_overlap(gene_lists)

  # Add a column to results with the group name
  res <- .assign_groups(res, gene_lists)

  # If the factor order is null, then use the order of the gene_lists names
  # The null_name and target_name should be the first two in the list.
  if (is.null(factor_order)) {
    factor_order <- names(gene_lists)
  }

  # Set order of groups for plotting based on the factor_order variable
  res$Group <- factor(res$Group, levels = factor_order)

  # Get the counts of genes in each gene lists
  gene_counts <- sapply(gene_lists, length)

  # Format the data better
  gene_counts <- sapply(gene_counts,
                        prettyNum,
                        big.mark = ",",
                        scientific = FALSE)

  # Create legend labels with counts
  legend_labels <- paste0(factor_order, ": ", gene_counts)

  # Use a stats test to determine significance
  test_results <- ecdf_stat_test(
    res,
    target_list = gene_lists[target_name],
    background_list = gene_lists[null_name],
    stats_test = stats_test,
    alternative = alternative
  )

  # Create ggplot
  ecdf_plot <- .plot_ecdf(
    res = res,
    palette = palette,
    title = title,
    legend_labels = legend_labels,
    test_results = test_results,
    theme = theme,
    x_lims = x_lims,
    output_filename = output_filename,
    height = height,
    width = width,
    dpi = dpi
  )

  return(list(plot = ecdf_plot,
              stats = test_results))
}

#' Check if input gene lists overlap
#'
#' @param gene_lists A list of gene lists.
#' example: list(c("gene1", "gene2"), c("gene1"))
#'
#' @return Warning if gene sets overlap.
#' @export
#'
#' @examples
#' # Overlap
#' SeedMatchR:::.check_gene_list_overlap(list(c("gene1", "gene2"), c("gene1")))
#'
#' #No overlap
#' SeedMatchR::.check_gene_list_overlap(list(c("gene1", "gene2"), c("gene3")))
.check_gene_list_overlap <- function(gene_lists) {
  if (any(duplicated(unlist(gene_lists)))) {
    warning("Gene lists overlap!")
  }
}

.assign_groups <- function(res, gene_lists) {
  res$Group  <- "NA"

  # For every group name in names(gene.list), label rows based on  list name
  for (group in names(gene_lists)) {
    res$Group <- ifelse(res$gene_id %in% unlist(gene_lists[group]),
                        group,
                        res$Group)
  }

  return(res)
}

#' Test for differences in log2(Fold Change) ECDFs between two gene lists
#'
#' @description This function uses the `stats` or `DTS` package to test the ECDF
#' of log2(Fold Changes) between two groups based on differential expression analysis.
#'
#' The inputs of this function are a DE results `data.frame` and two sets of
#' gene IDs called `target_list` and `background_list`. The functions will look
#' for a column called `log2FoldChange` in the data frame.
#'
#' The user can choose 4 tests for comparing the ECDFs of two gene sets.
#' Options: `KS`, `Wilcoxen`, `Wass`, or `DTS`.
#'
#' When looking at siRNA off-target effects use alternative = "greater".
#' If you are looking at effects on miRNA targets, use alternative = "lesser".
#' If you do not know which direction of change will occur, use alternative =
#' "two.sided".
#'
#' @param res Input results file data frame.
#' @param target_list Target set of genes.
#' @param background_list Set of genes that will be used as background
#' distribution.
#' @param stats_test Stats test to use. Options: KS, Wilcoxen, Wass, and DTS.
#' @param alternative The null hypothesis to test. Options: greater, less,
#' two.sided.
#'
#' @return A vector containing the test statistic and p-value.
#' @export
#' @examplesIf interactive()
#' library(dplyr)
#'
#' guide.seq = "UUAUAGAGCAAGAACACUGUUUU"
#'
#' anno.db = load_annotations("rnor7")
#'
#' # Load test data
#' get_example_data("sirna")
#'
#' sirna.data = load_example_data("sirna")
#'
#' res <- sirna.data$Schlegel_2022_Ttr_D1_30mkg
#'
#' # Filter DESeq2 results for SeedMatchR
#' res = filter_res(res, fdr.cutoff=1, fc.cutoff=0, rm.na.log2fc = TRUE)
#'
#' res = SeedMatchR(res = res, gtf = anno.db$gtf, seqs = anno.db$seqs,
#' sequence = guide.seq, seed.name = "mer7m8", tx.id.col= FALSE)
#'
#' # Gene set 1
#' mer7m8.list = res$gene_id[res$mer7m8 >= 1]
#'
#' # Gene set 2
#' background.list = res$gene_id[!(res$mer7m8 %in% mer7m8.list)]
#'
#' ecdf.res = SeedMatchR::ecdf_stat_test(res, mer7m8.list, background.list)
ecdf_stat_test <- function(res,
                           target_list,
                           background_list,
                           stats_test = c("KS", "Wilcoxen", "Wass", "DTS"),
                           alternative = c("greater", "two.sided", "less")) {
  alternative <- match.arg(alternative)
  stats_test <- match.arg(stats_test)

  # Get fold changes for the target group
  target_fc <-
    res$log2FoldChange[res$gene_id %in% unlist(target_list)]

  # Get fold chagnes for the background group
  background_fc <-
    res$log2FoldChange[res$gene_id %in% unlist(background_list)]

  if (stats_test == "KS") {
    scored <- stats::ks.test(target_fc,
                             background_fc,
                             alternative = alternative)

  } else if (stats_test == "DTS") {
    scored <- twosamples::dts_test(target_fc, background_fc)

    # These entries are added to the list with matching names
    scored$statistic <- scored[[1]]
    scored$p.value <- scored[[2]]

  } else if (stats_test == "Wass") {
    scored <- wass_dist(target_fc, background_fc)

    # These entries are added to the list with matching names
    scored$p.value <- NULL
    scored$statistic <- scored$neg.cens.auc

  } else {
    scored <- stats::wilcox.test(target_fc,
                                 background_fc,
                                 alternative = alternative)
  }

  return(scored)
}

#' Use ggplot2 to plot the ecdf for log2 fold changes
#'
#' @param res DE results data.frame.
#' @param palette Palette to use for line colors.
#' @param title Plot tile.
#' @param legend_labels Labels for legend entries.
#' @param test_results Results from the ecdf_stat_test function.
#' @param theme The ggplot2 theme parameters.
#' @param x_lims The x-axis limits for plotting.
#' @param output_filename Optional output filename.
#' @param height Height of output image in inches.
#' @param width Width of output image in inches.
#' @param dpi DPI of output image.
#'
#' @return A ggplot of the log2(fold changes) ECDF.
#' @export
#'
#' @examplesIf interactive()
#'
.plot_ecdf <- function(res,
                       palette = SeedMatchR.palette,
                       title,
                       legend_labels,
                       test_results,
                       theme = SeedMatchR.theme,
                       x_lims = c(-1, 1),
                       output_filename = NULL,
                       height = 5,
                       width = 5,
                       dpi = 300) {
  if (is.null(test_results$p.value)) {
    pvalue_label <- "Undefined"

  } else if (test_results$p.value < .01 &&
             test_results$p.value >= .001) {
    pvalue_label <- round(test_results$p.value, 3)

  } else if (test_results$p.value > .01) {
    pvalue_label <- round(test_results$p.value, 2)

  } else if (test_results$p.value < .001) {
    pvalue_label <- "p<0.001"
  }

  out_plot <- ggplot2::ggplot(
    as.data.frame(res),
    ggplot2::aes(
      x = log2FoldChange,
      group = Group,
      col = Group,
      geom = "step"
    )
  ) +
    ggplot2::stat_ecdf(linewidth = 1.5) +
    ggplot2::scale_colour_manual(values = palette, labels = legend_labels) +
    ggplot2::labs(
      title = title,
      x = "Log2(Fold Change)",
      y = "Cumulative Fraction",
      subtitle = paste0(
        "Statistic: ",
        round(test_results$statistic, 3),
        "\n",
        "p-value: ",
        pvalue_label
      )
    ) +
    SeedMatchR.theme +
    ggplot2::coord_cartesian(xlim = x_lims)

  # If there is an output name provided, save the figure
  if (!is.null(output_filename)) {
    ggplot2::ggsave(
      output_filename,
      out_plot,
      device = "pdf",
      height = height,
      width = width,
      units = "in",
      dpi = dpi
    )
  }

  return(out_plot)
}
