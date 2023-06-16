#' Check if input gene lists overlap
#'
#' @param gene.lists A list of gene lists.
#' example: list(c("gene1", "gene2"), c("gene1"))
#'
#' @return Warning if gene sets overlap
#' @export
#'
#' @examples
#' # Overlap
#' check_gene_list_overlap(list(c("gene1", "gene2"), c("gene1")))
#'
#' #No overlap
#' check_gene_list_overlap(list(c("gene1", "gene2"), c("gene3")))
check_gene_list_overlap <- function(gene.lists){
  if (any(duplicated(unlist(gene.lists)))){
    warning("Gene lists overlap!")
  }
}
