utils::globalVariables(c("baseMean", "padj","log2FoldChange", "pvalue"))

#' Filter differential expression results
#'
#' @description Filter differential expression results for use with search
#' and counting functions.
#'
#' The filtering criteria are:
#'
#' Filter out genes that are not expressed or counted at all:
#' baseMean == NA | pvalue == NA | log2FoldChange == NA
#'
#' Filter out genes with extreme outliers: pvalue = NA and padj = NA
#'
#' Filter out genes that have been excluded by independent filtering. padj = NA
#'
#' Filter results by the `fdr_cutoff`
#'
#' Filter the results by the `log2FoldChange`
#'
#' Filter the results by the `baseMean`
#'
#' @param res The differential expression results as a data frame
#' @param fdr_cutoff The false discovery rate cutoff to use.
#' @param fc_cutoff The fold change cutoff to use. The absolute
#' value will be used as the cutoff and values greater-than-or-equal-to will
#' be kept.
#' @param basemean_cutoff The minimum baseMean expression cutoff
#'
#' @return A filtered results data frame
#'
#' @export
#' @examplesIf interactive()
#' # Load test data
#' get_example_data("sirna")
#'
#' sirna.data = load_example_data("sirna")
#'
#' res <- sirna.data$Schlegel_2022_Ttr_D1_30mkg
#'
#' # Filter differential expression results for SeedMatchR
#' res = filter_res(res, fdr_cutoff=1, fc_cutoff=0)
filter_res <- function(res,
                         fdr_cutoff = 1,
                         fc_cutoff = 0,
                         basemean_cutoff = 30) {
  stopifnot("gene_id" %in% colnames(res),
            "log2FoldChange" %in% colnames(res),
            "padj" %in% colnames(res),
            "pvalue" %in% colnames(res),
            "baseMean" %in% colnames(res))

  # Filter out genes that are not expressed or counted at all: baseMean = 0,
  # pvalue = NA,  log2FoldChange = NA
  res <- res %>% dplyr::filter(!(is.na(padj) | is.na(pvalue) | is.na(log2FoldChange)))

  # Filter out genes with extreme outliers: pvalue = NA and padj = NA
  res <- res %>% dplyr::filter(!(is.na(pvalue) | is.na(padj)))

  # Filter out genes that have been excluded by independent filtering
  res <- res %>% dplyr::filter(!(is.na(padj)))

  # Filter results by the fdr_cutoff
  res <- res %>% dplyr::filter(padj <= fdr_cutoff)

  # Filter the results by the log2FoldChange
  res <- res %>% dplyr::filter(abs(log2FoldChange) >= fc_cutoff)

  # Filter the results by the baseMean
  res <- res %>% dplyr::filter(baseMean >= basemean_cutoff)

  return(res)

}
