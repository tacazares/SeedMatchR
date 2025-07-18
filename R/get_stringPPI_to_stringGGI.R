utils::globalVariables(c("combined_score"))

#' @title get_stringPPI_to_stringGGI
#'
#' @description This function retrieves the protein-protein interaction network from STRINGdb and maps the protein ids to gene ids using Biomart.
#'
#' @param sp_number The species number for the STRINGdb object
#' @param p The percentile threshold for the STRINGdb network
#' @param sp The species for the Biomart object
#'
#' @return A protein-protein interaction network with gene ids
#' @export
#' @examplesIf interactive()
get_stringPPI_to_stringGGI = function(sp_number, p = 0.95, sp = 'rat'){
  string_db <- STRINGdb::STRINGdb$new(species=sp_number)
  human_graph <- string_db$get_graph()

  edge.scores <- igraph::E(human_graph)$combined_score
  ninetyth.percentile <- stats::quantile(edge.scores, p)
  thresh <- data.frame(name= paste0(p*10,'th percentile'),
                       val=ninetyth.percentile)
  human_graph <- igraph::subgraph.edges(human_graph,
                                        igraph::E(human_graph)[combined_score > ninetyth.percentile])

  # 3. create adjacency matrix
  adj_matrix <- igraph::as_adjacency_matrix(human_graph)


  # 4. map gene ids to protein ids

  ### get gene/protein ids via Biomart
  if (sp == 'rat'){
    mart=biomaRt::useMart("ensembl",
                 dataset='rnorvegicus_gene_ensembl')
  }else if (sp == 'human'){
    mart=biomaRt::useMart(host = 'grch37.ensembl.org',
                 biomart='ENSEMBL_MART_ENSEMBL',
                 dataset='hsapiens_gene_ensembl')
  }


  ### extract protein ids from the human network
  protein_ids <- sapply(strsplit(rownames(adj_matrix), '\\.'),
                        function(x) x[2])

  ### get protein to gene id mappings
  mart_results <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                       "ensembl_peptide_id"),
                        filters = "ensembl_peptide_id", values = protein_ids,
                        mart = mart)

  ### replace protein ids with gene ids
  ix <- match(protein_ids, mart_results$ensembl_peptide_id)
  ix <- ix[!is.na(ix)]

  newnames <- protein_ids
  newnames[match(mart_results[ix,'ensembl_peptide_id'], newnames)] <-
    mart_results[ix, 'ensembl_gene_id']
  rownames(adj_matrix) <- newnames
  colnames(adj_matrix) <- newnames

  ppi <- adj_matrix[!duplicated(newnames), !duplicated(newnames)]
  nullrows <- Matrix::rowSums(ppi)==0
  ppi <- ppi[!nullrows,!nullrows] ## ppi is the network with gene ids

  return(ppi)
}
