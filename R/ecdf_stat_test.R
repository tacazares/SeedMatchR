#' Test for differences in log2(Fold Change) ECDFs between two gene lists using
#' the `twosamples` package
#'
#' @description This function uses the `twosamples` package to test the ECDF
#' of log2(Fold Changes) between two groups based on DESeq2 analysis.
#'
#' The inputs of this function are a DESeq2 results `data.frame` and two sets of
#' gene IDs called `gene.list1` and `gene.list2`. The functions will look for a
#' column called `log2FoldChange` in the dataframe.
#'
#' The user can choose 1 of 6 tests for comparing the ECDFs of two gene sets.
#' Options: `KS`, `Kuiper`, `DTS`, `CVM`, `AD`, `Wass`
#'
#' @param res Input results file data frame
#' @param gene.list1 Gene list 1: Usually null distribution
#' @param gene.list2 Gene list 2: Target set of genes
#' @param stats.test Stats test to use. Options: KS, Kuiper, DTS, CVM, AD, Wass
#' @param nboots Number of bootstraps to use for generating a null
#' @param keep.boots Whether to store the boostrap values
#'
#' @return A vector containing the dstat and pvalue
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
ecdf_stat_test <- function(res,
                           gene.list1,
                           gene.list2,
                           stats.test="KS",
                           nboots=10000,
                           keep.boots = TRUE){
  a = res$log2FoldChange[res$gene_id %in% unlist(gene.list1)]  #Null group

  b = res$log2FoldChange[res$gene_id %in% unlist(gene.list2)]  # A target

  # The test returns a class object that is not correctly interpreted by an
  # ifelse() function call!
  if(stats.test == "KS"){
    scored = twosamples::ks_test(a, b, nboots=nboots,
                                 keep.boots = keep.boots)
  } else if (stats.test == "Kuiper") {
    scored = twosamples::kuiper_test(a, b, nboots=nboots,
                                     keep.boots = keep.boots)
  } else if (stats.test == "DTS") {
    scored = twosamples::dts_test(a, b, nboots=nboots,
                                  keep.boots = keep.boots)
  } else if (stats.test == "CVM") {
    scored = twosamples::cvm_test(a, b, nboots=nboots,
                                  keep.boots = keep.boots)
  } else if (stats.test == "AD") {
    scored = twosamples::ad_test(a, b, nboots=nboots,
                                 keep.boots = keep.boots)
  } else if (stats.test == "Wass") {
    scored = twosamples::wass_test(a, b, nboots=nboots,
                                   keep.boots = keep.boots)
  } else {
    stop("Incorrect input!")
  }

  return(scored)
}
