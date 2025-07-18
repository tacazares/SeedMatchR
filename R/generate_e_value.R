utils::globalVariables(c("random_siRNA_results_overall", "symbol"))

#' generate_e_value
#'
#' @description This function calculates the e-value for each gene in the result_output dataframe. The e-value is calculated as the number of times the random siRNA sequences have a higher number of matches than the actual siRNA sequences.
#'
#' @param result_output A dataframe of siRNA sequences and their respective number of matches
#' @param random_seq_results A dataframe of random siRNA sequences and their respective number of matches
#'
#' @return A dataframe with the e-value for each gene
#' @export
generate_e_value = function(result_output, random_seq_results = random_siRNA_results_overall){

  temp_df = result_output %>% dplyr::filter(symbol %in% random_seq_results$symbol)
  # adding log transformation to normalize distribution
  temp_random_matrix = unname(log(as.matrix(random_seq_results[,2:ncol(random_seq_results)])))
  temp_df$e_value = 0
  for(k in temp_df$symbol){

    temp_df$e_value[which(temp_df$symbol == k)] = length(which(temp_random_matrix[which(random_seq_results$symbol == k),] >= log(temp_df$overall_matches[which(temp_df$symbol == k)])))/ncol(temp_random_matrix)

  }

  return(temp_df)
}
