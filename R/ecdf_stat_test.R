#' Test for differences in log2(Fold Change) ECDFs between two gene lists using
#' the `stats` package
#'
#' @description This function uses the `stats` package to test the ECDF
#' of log2(Fold Changes) between two groups based on DESeq2 analysis.
#'
#' The inputs of this function are a DESeq2 results `data.frame` and two sets of
#' gene IDs called `gene.list1` and `gene.list2`. The functions will look for a
#' column called `log2FoldChange` in the dataframe.
#'
#' @param res Input results file data frame
#' @param gene.list1 Gene list 1: Usually null distribution
#' @param gene.list2 Gene list 2: Target set of genes
#' @param stats.test Stats test to use. Options: KS or Wilcoxen
#' @param alternative The alternative hypothesis to test. Options: greater, less, two.sided
#'
#' @return A vector containing the dstat and pvalue
#'
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
#' sirna.data = load_example_data("sirna")
#'
#' res <- sirna.data$Schlegel_2022_Ttr_D1_30mkg
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
                           alternative="greater"){
  a = res$log2FoldChange[res$gene_id %in% unlist(gene.list1)] # Target group

  b = res$log2FoldChange[res$gene_id %in% unlist(gene.list2)] # Null group

  if(stats.test == "KS"){
    scored = stats::ks.test(a, b, alternative = alternative)

    } else if (stats.test == "Wilcoxen"){

      scored = stats::wilcox.test(a, b, alternative = alternative)

      } else {

        stop("Incorrect input!")
  }

  return(scored)
}
