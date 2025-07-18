utils::globalVariables(c("padj", "gene_id","log2FoldChange"))


#' Benchmark seed match predictions
#'
#' @param res A `data.frame` of `DESeq2` results.
#' @param preds_col Column name containing seed match counts.
#' @param padj Adjusted p-value cutoff for gold standard.
#' @param l2fc Log2(fold change) cutoff for gold standard.
#'
#' @return Two figures:
#'* A precision recall curve
#'* A receiver operator curve
#' @export
#'
#' @examplesIf interactive()
#' benchmark_preds(res)
benchmark_preds <- function(res,
                            preds_col = "mer7m8",
                            padj = .15,
                            l2fc = 0) {
  de_genes <- res %>%
    dplyr::filter(padj <= padj, log2FoldChange < l2fc) %>%
    dplyr::select(gene_id) %>%
    dplyr::distinct() %>%
    dplyr::pull()

  res$GS.label <- ifelse(res$gene_id %in% de_genes, 1, 0)

  preds <- res %>% dplyr::select(preds_col) %>% dplyr::pull()

  pr <- PRROC::pr.curve(scores.class0 = preds,
                        weights.class0 = res$GS.label,
                        curve = TRUE)

  roc <- PRROC::roc.curve(scores.class0 = preds,
                          scores.class1 = res$GS.label,
                          curve = TRUE)

  return(list("PRC" = pr,
              "ROC" = roc))
}
