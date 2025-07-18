#' Create a graph of protein-protein interactions from stringDB
#'
#' @description This function creates a graph object from the protein-protein interaction network and the result_output dataframe.
#'
#' @param string_ppi The protein-protein interaction network.
#' @param res The result_output data frame.
#'
#' @return A graph object with the feature attribute.
#' @export
create_graph <- function(string_ppi, res){
  g = igraph::graph_from_adjacency_matrix(string_ppi)

  igraph::V(g)$feature <- 0

  f = igraph::V(g)$feature

  names(f) = names(igraph::V(g))

  for( gene in names(f)){
    if(gene %in% res$gene_id){
      f[gene] = res$overall_matches[which(res$gene_id == gene)]
    }
  }

  igraph::V(g)$feature = f

  return(g)
}
