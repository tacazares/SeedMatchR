utils::globalVariables(c("mismatch_", "limit","TX_ID","TranscriptLocation","CompensatoryLength","check_d_tl","transcription_network", "Compensatory_Start"))

relative_position_seq = function(compensatory_region, guide.seq){
  x = as.data.frame(Biostrings::matchPattern(compensatory_region,gsub('U','T',guide.seq)))
  x_new = list(start = paste0(x$start,collapse = ','),
               end = paste0(x$end,collapse = ','))
  return(x_new)
}

#' @title get_compensatory_sequences
#' @name get_compensatory_sequences
#' @description Function to get 3' compensatory sequences for a given siRNA sequence.
#' @param guide.seq A character string representing the guide sequence.
#' @param seed A character string representing the seed sequence.
#' @param sequences A DNAStringSet object containing the sequences of the transcripts.
#' @param compensatory_region_length_min An integer representing the minimum length of the compensatory region.
#' @param compensatory_region_length_max An integer representing the maximum length of the compensatory region.
#' @param limit An integer representing the limit for sequence update.
#' @param mismatch_ An integer representing the maximum number of mismatches allowed in the seed match.
#' @param mismatch_compensatory_region_ An integer representing the maximum number of mismatches allowed in the compensatory region.
#' @return A data frame containing the compensatory regions.
#' @examplesIf interactive()
#'
#' @export
#'
get_compensatory_sequences = function(guide.seq,
                                      seed,
                                      sequences,
                                      compensatory_region_length_min = 6,
                                      compensatory_region_length_max = NA,
                                      limit=100,
                                      mismatch_ = 1,
                                      mismatch_compensatory_region_ = 0) {
  seed = gsub("U", "T", seed)

  non_can_seed_match_transcripts <- names(sequences)[which(Biostrings::vcountPattern(seed,
                                                                                     sequences,
                                                                                     max.mismatch = mismatch_,
                                                                                     with.indel = FALSE) > 0)]

  non_can_seed_match_transcripts <- sequences[non_can_seed_match_transcripts]

  stopifnot(is.character(guide.seq), nchar(guide.seq) >= 8)
  comp_seq <- substr(guide.seq, 8, nchar(guide.seq))

  if(is.na(compensatory_region_length_max)){
    compensatory_region_length_max <- nchar(comp_seq)
  }

  # All possible substrings of indexes of the compensatory region
  comp_seqs <- c()
  for (comp_len in compensatory_region_length_min:compensatory_region_length_max){
    comp_seqs <- c(comp_seqs, sapply(1:(nchar(comp_seq) - comp_len + 1), function(i){ return(gsub(" ","",substr(comp_seq,i,i+comp_len-1)))}))
  }

  # Data Frame Initialization
  compensatory_regions_captured <- data.frame("Seed" = character(),
                                             "Compensatory_Region" = character(),
                                             "TX_ID" = character(),
                                             "Distance" = numeric(),
                                             "Seed_Region_in_Guide" = character(),
                                             "Seed_Start" = numeric(),
                                             "Seed_End" = numeric(),
                                             "Compensatory_Start" = numeric(),
                                             "Compensatory_End" = numeric(),
                                             "Compensatory_Region_in_Guide" = character())
  if (nrow(compensatory_regions_captured) == 0) {
  message("No compensatory regions found.")
  # Return an empty data frame with the expected columns
  expected_cols <- c("Seed_in_Transcript", "Compensatory_Region_in_Transcript", "TX_ID", "Distance",
                     "Seed_Region_in_Guide", "Seed_Start", "Seed_End", "Compensatory_Start",
                     "Compensatory_End", "Compensatory_Region_in_Guide", "mismatch_in_seed",
                     "GC_Content", "CompensatoryLength", "relative_position_guide_seq_start",
                     "relative_position_guide_seq_end")
  return(as.data.frame(matrix(ncol = length(expected_cols), nrow = 0, dimnames = list(NULL, expected_cols))))
}
  else{
  seed_matches = as.data.frame(unlist(Biostrings::vmatchPattern(as.character(seed), non_can_seed_match_transcripts, max.mismatch = mismatch_)))
  seed_matches <- stats::setNames(seed_matches, c("seed_start","seed_end","seed_width", "gene_name"))

  # LOOP OVER COMPENSATORY SEQS
  k <- 1
  for (seq in comp_seqs){

    rev_seq <- Biostrings::reverseComplement(Biostrings::DNAString(gsub("U", "T", seq)))
    tmp = Biostrings::vmatchPattern(as.character(rev_seq), non_can_seed_match_transcripts, max.mismatch = mismatch_compensatory_region_)  # subset of non_can_seed_matches matching kmer
    rev_seq_matches <- as.data.frame(unlist(tmp))
    rev_seq_matches <- stats::setNames(rev_seq_matches, c("start","end","width", "gene_name"))

    merged_matches <- merge(seed_matches, rev_seq_matches, by = "gene_name")
    # Reduced merges
    merged_matches <- merged_matches[ merged_matches$seed_start > merged_matches$end + 1 & merged_matches$seed_start < merged_matches$end + limit,]

    if (nrow(merged_matches) > 0){
      merged_matches$Seed_in_Transcript = as.character(seed)
      merged_matches$Compensatory_Region_in_Transcript = as.character(rev_seq)
      merged_matches$TX_ID = merged_matches$gene_name
      merged_matches$Distance = merged_matches$seed_start - merged_matches$end - 1
      merged_matches$Seed_Region_in_Guide = as.character(Biostrings::reverseComplement(seed))
      merged_matches$Seed_Start = merged_matches$seed_end
      merged_matches$Seed_End = merged_matches$seed_start
      merged_matches$Compensatory_Start = merged_matches$end
      merged_matches$Compensatory_End = merged_matches$start
      merged_matches$Compensatory_Region_in_Guide = as.character(Biostrings::reverseComplement(rev_seq))

      merged_matches$mismatch_in_seed = as.numeric(as.character(seed) != as.character(Biostrings::subseq(non_can_seed_match_transcripts[merged_matches$gene_name],
                                                                                                    start=merged_matches$seed_start,
                                                                                                    end=pmin(Biostrings::width(non_can_seed_match_transcripts[merged_matches$gene_name]), merged_matches$seed_end))))
      compensatory_regions_captured = rbind(compensatory_regions_captured, merged_matches)
    }

    k <- k + 1
  }

  # Clean up extra columns
  compensatory_regions_captured <- compensatory_regions_captured[c("Seed_in_Transcript", "Compensatory_Region_in_Transcript", "TX_ID", "Distance", "Seed_Region_in_Guide", "Seed_Start",
                                    "Seed_End", "Compensatory_Start", "Compensatory_End", "Compensatory_Region_in_Guide","mismatch_in_seed")]


    # DATA SUMMARY
    compensatory_regions_captured$GC_Content = sapply(compensatory_regions_captured$Compensatory_Region_in_Transcript, .gc_content)
    compensatory_regions_captured$CompensatoryLength = sapply(compensatory_regions_captured$Compensatory_Region_in_Transcript, nchar)
    compensatory_regions_captured = compensatory_regions_captured %>% dplyr::group_by(TX_ID, Compensatory_Start) %>%  dplyr::filter(CompensatoryLength == max(CompensatoryLength)) %>% dplyr::ungroup()
    compensatory_regions_captured$check_d_tl = compensatory_regions_captured$Compensatory_Start - compensatory_regions_captured$CompensatoryLength
    compensatory_regions_captured = compensatory_regions_captured %>% dplyr::group_by(TX_ID, check_d_tl) %>%  dplyr::filter(CompensatoryLength == max(CompensatoryLength)) %>% dplyr::ungroup()
    compensatory_regions_captured$check_d_tl = NULL
    compensatory_regions_captured$relative_position_guide_seq_start = sapply(compensatory_regions_captured$Compensatory_Region_in_Guide,function(x){return(relative_position_seq(x,guide.seq)$start)})
    compensatory_regions_captured$relative_position_guide_seq_end = sapply(compensatory_regions_captured$Compensatory_Region_in_Guide,function(x){return(relative_position_seq(x,guide.seq)$end)})
    return(compensatory_regions_captured)
  }


}


#' @title .update_sequence_at
#' @description Function to update the sequence library by extracting subsequences around the seed match.
#' @param sequences A DNAStringSet object containing the sequences of the transcripts.
#' @param seed A character string representing the seed sequence.
#' @param mismatch_allowed An integer representing the maximum number of mismatches allowed in the seed match.
#' @param limit_ An integer representing the limit for sequence update.
#' @return A DNAStringSet object containing the updated sequences and a data frame with the start sites of the seed matches.
.update_sequence_at <- function(sequences, seed, mismatch_allowed = mismatch_, limit_ = limit) {

    matches__ <- Biostrings::vmatchPattern(seed, sequences, max.mismatch = mismatch_allowed)
    matches__ = matches__[sapply(matches__, length) > 0]
    df_list_end_list <- lapply(matches__, function(x){return(max(IRanges::end(x)))})
    df_list_start_list <- lapply(matches__, function(x){return(min(IRanges::start(x)))})

    df_list_start = lapply(df_list_start_list,function(x){if(x>limit_){return(x-limit_)}else{return(1)}})
    df_list_end = as.integer(df_list_end_list)

    # Convert df_list_start and df_list to integer vectors
    df_list_start_int <- as.integer(df_list_start)
    df_list_int <- as.integer(df_list_start_int)
    seq_lengths <- Biostrings::width(sequences[names(df_list_start)])

    # Use mapply for vectorized operations
    combined_df <- mapply(function(seq, start, end, len) {
      return(Biostrings::subseq(seq, start, min(end, len)))
    }, sequences[names(df_list_start)], df_list_start_int, df_list_end, seq_lengths)
    combined_df = Biostrings::DNAStringSet(combined_df)
    start_sites = data.frame("transcript" = names(df_list_start), "SeedSite" = unname(as.numeric(df_list_start)))
    out = list(combined_df = combined_df, start_sites = start_sites)
    return(out)
}

#' @title .gc_content
#' @description Function to calculate the GC content of a given DNA sequence.
#' @param sequence A character string representing the DNA sequence.
#' @return A numeric value representing the GC content percentage of the sequence.
.gc_content <- function(sequence) {
  # Convert the sequence to uppercase to handle both cases
  stopifnot("Empty String provided" = length(sequence) > 0) |> try()

  sequence = toupper(sequence)

  # Count the number of G and C bases
  num_g = sum(strsplit(sequence, '')[[1]] == "G")
  num_c = sum(strsplit(sequence, '')[[1]] == "C")

  # Calculate the GC content
  gc_content = (num_g + num_c) / nchar(sequence) * 100

  return(gc_content)
}
