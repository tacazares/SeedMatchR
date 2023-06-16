#' Generate all possible sequence permutations given an alphabet and length
#'
#' @description This function is used to generate all possible sequences given
#' an alphabet (`alph`) and length of sequence (`n`).
#'
#' @param n number of characters in the sequence
#' @param alph the alphabet
#'
#' @return A data.frame containg all possible sequences. Each columns is a
#' position in the sequence.
#' @export
#' @examples
#' seq.df = generate_possible_seqs(n=8)
generate_possible_seqs <- function(n=8, alph = c("A", "U", "G", "C")){
  seed.set = expand.grid(replicate(n, alph, simplify = FALSE))

  return(seed.set)

}
