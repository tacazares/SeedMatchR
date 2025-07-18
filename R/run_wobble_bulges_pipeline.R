utils::globalVariables(c("bulges_lim", "nt", "bulge_allowance"))

#' Running wobble and bulge inclusive seedmatchR pipeline
#'
#' @description This function runs a wobble and bulge inclusive SeedMatchR pipeline.
#'
#' @param guide.seq A character string representing the guide sequence.
#' @param res A DESeq2 results data.frame.
#' @param features A DNAStringSet object with sequence information for features. The names of the sequences should be the transcript names.
#' @param bulges_lim An integer representing the number of target bulges allowed per sequence. Default is 4.
#' @param nt An integer representing the length from the start of siRNA where wobbles are not allowed. Default is 3.
#' @param silent A logical value to determine if the function should print progress. Default is FALSE.
#'
#' @return A data.frame with the overall matches for each seed.
#' @export
#' @examplesIf interactive()
#' result_df = run_wobble_bulges_pipeline(guide.seq = guide.seq, res = res, features = features)
run_wobble_bulges_pipeline = function(guide.seq, res = res, features = features, bulges_lim = 4, nt = 3, silent = FALSE){

  guide_seq_length <- nchar(guide.seq)

  combinations_w_tb = .get_all_combination_w_b(guide.seq[1], silent =FALSE)

  seed_seq_set_w_wobble = combinations_w_tb[[2]]
  seed_seq_set_w_wobble = sapply(seed_seq_set_w_wobble, function(x){return(gsub('T','U',as.character(x)))})

  all_seed_possible = sapply(seed_seq_set_w_wobble,function(x){return(.get_target_bulged_guides(x,bulges_lim, guide_seq_length))})
  all_seed_checked = sapply(all_seed_possible,function(x){return(check_bulges_nt(x,nt))})

  if(!silent){
    cat("\n Running for all seeds ....")
  }

  temp_out_mod <- NULL
  for(j in 1:length(all_seed_checked)){
    for(i in 1:length(all_seed_checked[[j]])){
      if(i == 1 || is.null(temp_out_mod)){
      #temp_out_mod = SeedMatchR(res, features, sequence = gsub("U", "T", all_seed_checked[[j]][i]), sirna.name = paste0("SeedMatches_",j),get_seed = TRUE,res.format = 'data.frame')
        temp_out_mod = data.frame("transcript" = names(features),"SeedMatches" = Biostrings::vcountPattern(pattern = Biostrings::DNAString(gsub("U", "T", all_seed_checked[[j]][i])),subject = features,fixed = FALSE))
        colnames(temp_out_mod) = c("transcript",paste0("SeedMatches_",j))
        }else{
        temp_out_mod[paste0("SeedMatches_",j)] = temp_out_mod[paste0("SeedMatches_",j)] + Biostrings::vcountPattern(pattern = Biostrings::DNAString(gsub("U", "T", all_seed_checked[[j]][i])),subject = features, fixed = FALSE)
      }
    }
    if(j > 1){
      out_mod_seed_checked = cbind(out_mod_seed_checked,temp_out_mod[paste0("SeedMatches_",j)])
    }else{
      out_mod_seed_checked = temp_out_mod
    }
  }
  if(length(grep('SeedMatches_',colnames(out_mod_seed_checked))) == 1){
    out_mod_seed_checked$overall_matches = out_mod_seed_checked[,grep('SeedMatches_',colnames(out_mod_seed_checked))]
  }else{
    out_mod_seed_checked$overall_matches = rowSums(out_mod_seed_checked[,grep('SeedMatches_',colnames(out_mod_seed_checked))])
  }

  if(!silent){
    cat("Completed \n")
  }

  return(out_mod_seed_checked)
}


#' @title get_all_combination_w_b
#'
#' @description This function generates all possible combinations of guide sequences with wobbles
#'
#' @param guide.seq The guide sequence
#' @param silent A boolean to suppress output
#'
#' @return A list of guide and seed sequences
#' @export
#' @examplesIf interactive()
.get_all_combination_w_b = function(guide.seq, silent = FALSE){

  split.string = grep('U',strsplit(guide.seq,'')[[1]])

  # Step 1 : Seed Scanning in UTR -------------------------------------------
  num_U = length(split.string)

  if (num_U > 3){
    total_combs =  factorial(num_U)/(factorial(4)*factorial(num_U-4))
    if(!silent){
      cat("\n Loading seeds with wobbles for ",total_combs," combinations \n")

    }
    input_list = split(as.data.frame(gtools::combinations(num_U,
                                                          4,
                                                          split.string)),
                       seq(nrow(gtools::combinations(num_U,
                                                     4,
                                                     split.string))))

  }else{
    if(!silent){

      cat("\n Less than three Uridine present \n")
    }
    total_combs = .find_factorial_internal(num_U = num_U)
    input_list = split(as.data.frame(gtools::combinations(num_U,num_U,grep('U',split.string)),seq(nrow(gtools::combinations(num_U,num_U, split.string)))))

  }

  guide_seq_set_w_wobble = sapply(input_list,
                                  function(x){
                                    temp_set = strsplit(guide.seq,'')[[1]]
                                    x = as.numeric(x)
                                    temp_set[x] = "Y"
                                    return(paste0(temp_set,collapse = ""))
                                  })

  seed.seq_set_w_wobble = sapply(guide_seq_set_w_wobble,function(x){return(SeedMatchR::get_seed(x,"mer7m8")$Target.Seq)})

  #Obtaining the unique seeds
  seed.seq_set_w_wobble_unique = c()

  for(i in 1:length(seed.seq_set_w_wobble)){
    if(length(seed.seq_set_w_wobble_unique) == 0 || !(as.character(seed.seq_set_w_wobble[i][[1]]) %in% seed.seq_set_w_wobble_unique)){
      seed.seq_set_w_wobble_unique = c(seed.seq_set_w_wobble_unique,as.character(seed.seq_set_w_wobble[i][[1]]))
    }
  }

  seed.seq_set_w_wobble_unique = sapply(seed.seq_set_w_wobble_unique,function(x){return(Biostrings::DNAString(x))})

  return(list("guide_combinations" = guide_seq_set_w_wobble,
              "seed_combinations" = seed.seq_set_w_wobble_unique))

}


#' @title find_factorial_internal
#'
#' @description This function calculates the sum of combinations for a given number of Uridines.
#'
#' @param num_U An integer representing the number of Uridines.
#'
#' @return The sum of combinations as an integer.
#' @export
#' @examplesIf interactive()
.find_factorial_internal = function(num_U){
  out = 0
  for(x in 1:num_U){
    out = out +  factorial(num_U)/(factorial(x)*factorial(num_U-x))
  }
  return(out)
}

#' @title get_target_bulged_guides
#'
#' @description This function generates a set of RNA sequences with bulges inserted at various positions.
#'
#' @param guide.seq A character vector representing the guide sequence.
#' @param guide_seq_length Integer, length of guide strand sequence
#' @param bulges_lim An integer representing the maximum number of bulges allowed.
#'
#' @return An RNAStringSet object containing the sequences with bulges.
#' @export
#' @examplesIf interactive()
.get_target_bulged_guides <- function(guide.seq, guide_seq_length, bulges_lim = bulges_lim){

  sequence = as.character(guide.seq)
  sequence = strsplit(sequence,'')[[1]]

  sequence_with_bulge_set = c()
  for(i in 1:bulges_lim){
    bulge_allowance <- i
    sequence_with_bulge_set = c(sequence_with_bulge_set,sapply(bulge_allowance:(length(sequence) - bulge_allowance - 1),function(x){return(.sequence_with_tb(x,sequence, guide_seq_length))}))
  }

  sequence_with_bulge_set = Biostrings::RNAStringSet(sequence_with_bulge_set)
  return(sequence_with_bulge_set)
}

#' @title .sequence_with_tb
#'
#' @description This nested function generates a sequence with a bulge inserted at the specified start position.
#' @param bulge_start An integer representing the start position of the bulge.
#' @param guide_seq_length Integer, length of guide strand sequence
#' @param sequence Input sequence
#'
#' @return A character string representing the sequence with the bulge.
#' @export
#' @examplesIf interactive()
.sequence_with_tb = function(bulge_start, sequence, guide_seq_length) {
  if (bulge_start <= guide_seq_length){
    k = 1
  }else{
    k = bulge_start - guide_seq_length
  }

  if(bulge_start +  bulge_allowance + 1 + guide_seq_length < length(sequence)){
    l = bulge_start +  bulge_allowance + 1 + guide_seq_length
  }else{
    l = length(sequence)
  }
  return(paste0(c(sequence[k:bulge_start],rep("N",bulge_allowance), sequence[(bulge_start + bulge_allowance + 1):l]), collapse = ""))
}

#' @title check_bulges_nt
#' @description This function checks for sequences with bulges and filters them based on the presence of 'N' in the first `nt` positions.
#' @param sequenceStringSet An RNAStringSet object containing the sequences to be checked.
#' @param nt An integer representing the length from the start of the sequence where bulges are not allowed.
#' @param return_allowed A boolean indicating whether to return the allowed sequences (TRUE) or the removed sequences (FALSE).
#' @return An RNAStringSet object containing either the allowed or removed sequences based on the `return_allowed` parameter.
#' @export
#' @examplesIf interactive()
check_bulges_nt = function(sequenceStringSet, nt = nt,return_allowed = TRUE){

  allowed_sequenceStringSet = c()
  removed_sequence = c()
  for(i in 1:length(sequenceStringSet)){
    t = strsplit(as.character(sequenceStringSet[[i]]),'')[[1]][1:nt]
    if('N' %in% t){
      removed_sequence = c(removed_sequence,i)
    }else{
      allowed_sequenceStringSet = c(allowed_sequenceStringSet,i)
    }
  }
  if(return_allowed){
    return(sequenceStringSet[allowed_sequenceStringSet])
  }else{
    return(sequenceStringSet[removed_sequence])
  }
}

