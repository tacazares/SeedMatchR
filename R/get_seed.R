#' Get the target seed sequence given a canonical seed name and input sequence
#'
#' @description Given a sequence greater than 8 bp oriented 5' -> 3' and a seed
#' definition, this function will return an object containing seed-specific sequence
#' information. Users can input a custom seed name, but must provide the start
#' position (`start.pos`) and stop position (`stop.pos`) that define the
#' range of the seed sequence.
#'
#' Built-in options: `mer8`,  `mer7A1`, `mer7m8`, `mer6`
#'
#' Note: The seed definitions `mer8` and `mer7A1` force a U at position g1.
#' This results in an A in the target sequence being searched.
#'
#' @param guide.seq A character string greater than 8 bp and oriented 5'-> 3'.
#' @param seed.name The seed name of interest. Options: mer8,  mer7A1, mer7m8,
#' mer6. If not in the default list, the start.pos and stop.pos arguments
#' will be used to define the seed.
#' @param start.pos The start position for a custom seed definition
#' @param stop.pos The stop position for a custom seed definition
#'
#' @return  An object with the entries:
#'
#' * `Guide`: Input guide sequence. Input is expected to be RNA.
#' * `Seed.Name`: The seed name.
#' * `Seed.Seq.RNA`: The seed sequence as a `RNAString`
#' * `Seed.Seq.DNA`: The seed sequence as a `DNAString`
#' * `Target.Seq`: The target DNA sequence based on the reverse complement of
#'    the seed as a `DNAString`
#'
#' @export
#'
#' @examples
#' # Example Ttr from Schlegel et al. 2022
#' guide.seq = "UUAUAGAGCAAGAACACUGUUUU"
#'
#' # Get seed match
#' seed.seq = get_seed(guide.seq, "mer7m8")
get_seed <- function(guide.seq, seed.name = "mer7m8", start.pos = 1, stop.pos = 8){
    seed.names.avail = c("mer7m8", "mer8", "mer6", "mer7A1")

    # Check parameters
    stopifnot(is.character(guide.seq),
            nchar(guide.seq) >= 8)

    # Test whether the seed.name is already defined. If not then use the input start and stop coordinates
    if (seed.name %in% seed.names.avail){
        # Get the substring of the sequence for the seed g2-g8
        base.seed = Biostrings::RNAString(substring(guide.seq, 2, 8))

        # Get the reverse complement of the sequence and convert to DNAString object
        seed.rev.comp.dna = Biostrings::reverseComplement(Biostrings::DNAString(base.seed))

        # Depending on seed.name, return the seed sequence: assume there will always
        # be an A at g1.
        seed.target.seq = ifelse(seed.name == "mer8",
                           paste0(seed.rev.comp.dna, "A"),
                           ifelse(seed.name == "mer7A1",
                                  paste0(substr(seed.rev.comp.dna, 2, 7), "A"),
                                  ifelse(seed.name == "mer7m8",
                                         paste0(seed.rev.comp.dna),
                                         ifelse(seed.name == "mer6",
                                                paste0(substr(seed.rev.comp.dna, 2, 7)),
                                                paste0("You broke SeedMatchR")))))
         # Convert to DNAString objects
         seed.target.seq = Biostrings::DNAString(seed.target.seq)

    } else {
         message(paste0("No seed definition selected. Setting name as: ",
                      seed.name,
                      ". Using g",
                      start.pos,
                      ":g",
                      stop.pos,
                      " as start:stop"))

         rna.string = Biostrings::RNAString(guide.seq)

         base.seed = substr(rna.string, start.pos, stop.pos)

         seed.seq.rna = Biostrings::RNAString(base.seed)

         seed.seq.dna = Biostrings::DNAString(seed.seq.rna)

         seed.target.seq = Biostrings::reverseComplement(seed.seq.dna)
    }

    return(c(Guide = Biostrings::RNAString(guide.seq),
             Seed.Name = seed.name,
             Seed.Seq.RNA = Biostrings::RNAString(Biostrings::reverseComplement(seed.target.seq)),
             Seed.Seq.DNA = Biostrings::reverseComplement(seed.target.seq),
             Target.Seq = seed.target.seq))
}
