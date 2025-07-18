#' Find inter match distances for each gene
#'
#' @param gr_list A GRangesList of matches from SeedMatchR
#' @param feature.pos Which position to use for central point. Options: "fivePrimeEnd", "threePrimeEnd", "center"
#'
#' @return A list of calculated distances for each match in the GRanges object grouped by gene name.
#' @export
#'
#' @examplesIf interactive()
#' match_distances(gr)
#'
match_distances <- function(gr_list, feature.pos = c("fivePrimeEnd", "threePrimeEnd", "center")) {
  feature.pos = match.arg(feature.pos)

  # Initialize a list to store the distances for each Grange
  dists_per_grange <- list()

  # Loop over all Granges in the list
  for (i in 1:length(gr_list)) {
    # Get the current Grange
    gr <- gr_list[[i]]

    feature_coordinates <- switch(
      feature.pos,
      "center" = .get_centers(gr),
      "fivePrimeEnd" = .get_fivePrimeEnd(gr),
      "threePrimeEnd" = .get_threePrimeEnd(gr))

    # Initialize a vector to store the distances for this Grange
    dists <- c()

    # Loop over all pairs of centers in the Grange
    for (j in 1:(length(feature_coordinates) - 1)) {
      for (k in (j + 1):length(feature_coordinates)) {
        # Calculate the distance between the two centers
        dist <- abs(feature_coordinates[j] - feature_coordinates[k])

        # Add the distance to the vector
        dists <- c(dists, dist)
      }
    }

    # Add the distances for this Grange to the list
    dists_per_grange[[names(gr_list)[i]]] <- dists
  }

  #names(dists_per_grange) = names(gr)

  # Return the distances per Grange
  return(dists_per_grange)
}


.get_centers <- function(gr){
  centers <- BiocGenerics::start(gr) + BiocGenerics::width(gr) / 2

  return(centers)
}


.get_fivePrimeEnd <- function(gr){
  fivePrimeEnds = BiocGenerics::end(gr)

  return(fivePrimeEnds)
}


.get_threePrimeEnd <- function(gr){
  fivePrimeEnds = BiocGenerics::start(gr)

  return(fivePrimeEnds)
}
