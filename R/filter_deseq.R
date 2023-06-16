#' Filter DESEQ2 Results for SeedMatchR
#'
#' @description Filter \code{DESeqDataSet} results for use with seed matching
#' and counting functions.
#'
#' The filtering criteria are:
#'
#' Filter out genes that are not expressed or counted at all:
#' baseMean = 0 & pvalue = NA & log2FoldChange = NA
#'
#' Filter out genes that are expressed, but there is not difference across
#' groups: log2FoldChange = 0
#'
#' Filter out genes with extreme outliers: pvalue = NA and padj = NA
#'
#' Filter out genes that have been excluded by independent filtering. padj = NA
#'
#' Filter results by the `fdr.cutoff`
#'
#' Filter the results by the `log2FoldChange`
#'
#' Filter the results by the `baseMean`
#'
#' Remove NA gene_ids and log2FoldChange values
#'
#' @param res The DESEQ2 results as a data frame
#' @param fdr.cutoff The false discovery rate cutoff to use.
#' @param fc.cutoff The fold change cutoff to use. The absolute
#' value will be used as the cutoff and values greater-than-or-equal-to will
#' be kept.
#' @param rm.na.name Remove na values from the gene_name column
#' @param rm.na.log2fc Remove na values from the log2FoldChange column
#' @param baseMean.cutoff The minimum baseMean expression cutoff
#'
#' @return A modified DESEQ2 results table that has been filtered
#'
#' @export
#' @examplesIf interactive()
#' # Load test data
#' get_example_data("sirna")
#'
#' load_example_data("sirna")
#'
#' res <- Schlegel_2022_Ttr_D1_30mkg
#'
#' # Filter DESeq2 results for SeedMatchR
#' res = filter_deseq(res, fdr.cutoff=1, fc.cutoff=0, rm.na.log2fc = TRUE)
filter_deseq <- function(res,
                         fdr.cutoff=1,
                         fc.cutoff=0,
                         rm.na.name=FALSE,
                         rm.na.log2fc = FALSE,
                         baseMean.cutoff=0){
  # Filter out genes that are not expressed or counted at all: baseMean = 0 &
  # pvalue = NA & log2FoldChange = NA
  res <- res[!(res$baseMean == 0 & is.na(res$padj) & is.na(res$pvalue) &
                 is.na(res$log2FoldChange)),]

  # Filter out genes that are expressed, but there is not difference across
  # groups: log2FoldChange = 0
  res <- res[!res$log2FoldChange == 0,]

  # Filter out genes with extreme outliers: pvalue = NA and padj = NA
  res <- res[!(is.na(res$pvalue) & is.na(res$padj)),]

  # Filter out genes that have been excluded by independent filtering
  res <- res[!is.na(res$padj),]

  # Filter results by the fdr.cutoff
  res <- res[res$padj <= fdr.cutoff,]

  # Filter the results by the log2FoldChange
  res <- res[abs(res$log2FoldChange) >= fc.cutoff,]

  # Filter the results by the baseMean
  res <- res[res$baseMean >= baseMean.cutoff,]

  if (rm.na.name == TRUE){
    # Remove rows with no gene name
    res = res[!res$gene_name == "",]

    # Remove rows with NA in the gene name column
    res = res[!is.na(res$gene_name),]
  }

  if (rm.na.log2fc == TRUE){
    # Remove rows with NA values in the FC column
    res = res[!is.na(res$log2FoldChange),]
  }

  return(res)

}
