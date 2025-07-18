expected.names = c("Full-D0", "Full-D1", "Full-D2", "Full-D3",
                   "Full-D4", '18mer-D0', "18mer-D1", "18mer-D2",
                   "18mer-D3", "18mer-D4", "15mer-D4", "15mer-D3",
                   "15mer-D2", "15mer-D1",  "15mer-D0", "8mer", "7mer-m8",
                   "7mer-A1","6mer")

utils::globalVariables(c("Distance", "best.match", "Group", "gene_id",
                         "IN.DESEQ2", "log2FoldChange", "padj", expected.names))

#' Generate statistics based on sequence matches and DE data
#'
#' @param res A `data.frame` of `DESeq2` results
#' @param seqs A `Biostrings::DNAStringSet` of sequences to search
#' @param guide.seq Input RNA sequence as a character vector
#' @param group.name The group name to use for labeling rows
#' @param indel.bool Boolean in long format (TRUE, FALSE) for whether to include indels.
#' @param allow_wobbles Whether to allow G:U wobbles in the search.
#' @param log2fc.cut Log2(FC) threshold for considering a gene an off-target.
#' @param padj.cut  Adjust p-value threshold
#' @param fixed Whether IUPAC values need to match for ambiguous chracters.
#' @param seed.edit.distances An integer value indicating allowed edit distance for seed searches.
#'
#' @return A list containing:
#'            report: a `data.frame` of results
#'            table: a `gt::gt` table of results
#'            matches: A data.frame of matches
#'
#' @export
#'
#' @examplesIf interactive()
#' #' guide.seq = "UUAUAGAGCAAGAACACUGUUUU"
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
#' SeedMatchReport(res, anno.db$seqs, anno.db$gtf, guide.seq, "si-Ttr")
SeedMatchReport <- function(res,
                            seqs,
                            guide.seq,
                            group.name = "SeedMatchReport",
                            indel.bool = TRUE,
                            allow_wobbles = FALSE,
                            log2fc.cut = -0.7655347,
                            padj.cut = 0.05,
                            seed.edit.distances = 0,
                            fixed = T){

  # Perform a full search of the annotation database across multiple patterns
  match_matrix = full_search(guide.seq,
                             seqs,
                             group.name,
                             indel.bool = indel.bool,
                             allow_wobbles = allow_wobbles,
                             seed.edit.distances = seed.edit.distances,
                             fixed = fixed)

  # add column to matches based on the sequence names
  match_matrix$gene_id = names(seqs)

  # Keep only sequences in the input sequence file
  res = res %>% dplyr::filter(gene_id %in% names(seqs))

  # Create a column that labels genes as detected or undetected
  match_matrix$IN.DESEQ2 = ifelse(match_matrix$gene_id %in% res$gene_id,
                                  "Detected",
                                  "Undetected")

  # Map log2FoldChange, padj, and baseMean from DE data frame
  match_matrix$log2FoldChange = res$log2FoldChange[match(match_matrix$gene_id,
                                                         res$gene_id)]

  match_matrix$padj = res$padj[match(match_matrix$gene_id, res$gene_id)]

  match_matrix$baseMean =  res$baseMean[match(match_matrix$gene_id, res$gene_id)]

  in.silico.count = match_matrix %>%
                    dplyr::group_by(best.match) %>%
                    dplyr::summarize("In silico predictions" = dplyr::n())

  detected.counts = match_matrix %>%
                    dplyr::filter(IN.DESEQ2 == "Detected") %>%
                    dplyr::group_by(best.match) %>%
                    dplyr::summarize("Expressed predictions" = dplyr::n())

  DE.counts = match_matrix %>%
              dplyr::filter(log2FoldChange <= log2fc.cut,
                            padj <= padj.cut,
                            best.match != "No.matches") %>%
              dplyr::group_by(best.match) %>%
              dplyr::summarize("Off-target predictions" = dplyr::n())

  data_list = list(in.silico.count, detected.counts, DE.counts)

  merged.counts <- Reduce(function(x, y) merge(x, y, by = "best.match", all = TRUE), data_list)

  merged.counts = merged.counts %>%
    tibble::column_to_rownames("best.match")

  merged.counts = merged.counts[rownames(merged.counts) != "No.matches",]

  missing.rows = expected.names[!(expected.names %in% rownames(merged.counts))]

  for (miss.name in missing.rows){
    merged.counts = rbind(merged.counts, miss.name = c(0,0,0))

    rownames(merged.counts)[nrow(merged.counts)] <- miss.name
  }

  merged.counts[is.na(merged.counts)] <- 0

  merged.counts = rbind(merged.counts, "Total" = colSums(merged.counts))

  merged.counts$`% off-target` = (merged.counts$`Off-target predictions` /
                                       merged.counts$`Expressed predictions`) * 100

  merged.counts = merged.counts[c(expected.names, "Total"),]

  merged.counts[is.na(merged.counts)] <- 0

  res.df = as.data.frame(t(merged.counts))

  match.report = tibble::rownames_to_column(res.df, "Category")

  match.report$Group = group.name

  match.table = .make_table(match.report)

  return(list(report = match.report,
              table = match.table,
              matches = match_matrix))
}

.make_table <- function(match.report){
  match.table = match.report |>
    dplyr::group_by(Group) |>
    gt::gt(rowname_col = "Category", id = "mygt") |>
    gt::fmt_number(columns = c(expected.names, "Total"),
                   rows = c("In silico predictions",
                            "Expressed predictions",
                            "Off-target predictions"),
                   decimals = 0)  |>
    gt::fmt_percent(columns = c(expected.names, "Total"),
                    rows = c("% off-target"),
                    decimals = 2,
                    scale_values = F) |>
    gt::tab_header(title = "In-silico siRNA Binding Prediction",
                   subtitle = "Identifying siRNA hits in the transcriptome") |>
    gt::tab_spanner(
      label = "Full Guide Strand (g2:g23)",
      columns = c(`Full-D0`, `Full-D1`, `Full-D2`, `Full-D3`, `Full-D4`)
    ) |>
    gt::tab_spanner(
      label = "18-mer (g2:g19)",
      columns = c(`18mer-D0`, `18mer-D1`, `18mer-D2`, `18mer-D3`, `18mer-D4`)
    ) |>
    gt::tab_spanner(
      label = "15-mer (g2:g19)",
      columns = c(`15mer-D0`, `15mer-D1`, `15mer-D2`, `15mer-D3`, `15mer-D4`)
    ) |>
    gt::data_color(
      columns = expected.names,
      rows = gt::contains(c("% off-target")),
      method = "numeric",
      palette = c("white", "red"),
      domain = c(0, 100)
    ) |>
    gt::cols_label(
      `Full-D0` = "D0",
      `Full-D1` = "D1",
      `Full-D2` = "D2",
      `Full-D3` = "D3",
      `Full-D4` = "D4",
      `18mer-D0` = "D0",
      `18mer-D1` = "D1",
      `18mer-D2` = "D2",
      `18mer-D3` = "D3",
      `18mer-D4` = "D4",
      `15mer-D0` = "D0",
      `15mer-D1` = "D1",
      `15mer-D2` = "D2",
      `15mer-D3` = "D3",
      `15mer-D4` = "D4",
    ) |>
    gt::opt_stylize(style = 3, add_row_striping = T)

  return(match.table)
}
