#' Get the target sequence for a given siRNA input sequence
#'
#' @description Given a sequence greater than 8 bp oriented 5' -> 3' and a seed
#' definition (or custom definition), this function will return an object
#' containing target sequence information.
#'
#' Built-in options: `mer8`,  `mer7A1`, `mer7m8`, `mer6`, `Full`, and `18mer`
#'
#' Custom Definitions: Users can input a custom seed name, but must provide the
#' start position (`start.pos`) and stop position (`stop.pos`) that define the
#' range of the seed sequence.
#'
#' Note: `mer8` and `mer7A1` force a U at position g1. This results in an A in the target
#' sequence being searched. This is based on historical data and observations
#' from miRNA biology. If you do not wish to enforce this, provide the start,
#' stop, and custom name to use the original sequence.
#'
#' @param guide.seq A character string greater than 8 bp and oriented 5'-> 3'.
#' @param seed.name The seed name of interest. Options: mer8,  mer7A1, mer7m8,
#' mer6, 18-mer,  15-mer, or Full. If not in the default list, the start.pos and stop.pos arguments
#' will be used to define the seed.
#' @param start.pos The start position for a custom seed definition. 1-based
#' @param stop.pos The stop position for a custom seed definition.
#' @param allow_wobbles If True, allow G:U wobbles by replacing U with Y.
#'
#' @return  An object with the entries:
#'
#' * `Guide`: Input guide sequence. Input is expected to be RNA.
#' * `Seed.Name`: The seed name.
#' * `Seed.Seq.RNA`: The seed sequence as a `RNAString`
#' * `Seed.Seq.DNA`: The seed sequence as a `DNAString`
#' * `Target.Seq`: The target DNA sequence based on the reverse complement of
#'    the `Seed.Seq.DNA`
#'
#' @export
#'
#' @examples
#' # Example Ttr from Schlegel et al. 2022
#' guide.seq = "UUAUAGAGCAAGAACACUGUUUU"
#'
#' # Get seed match
#' seed.seq = get_seed(guide.seq = guide.seq, seed.name = "mer7m8")
get_seed <- function(guide.seq,
                     seed.name = "mer7m8",
                     start.pos = NULL,
                     stop.pos = NULL,
                     allow_wobbles = FALSE){
    # Check parameters
    stopifnot(is.character(guide.seq))

    cannonical.seeds = c("mer7m8", "mer8", "mer6", "mer7A1", "Full", "18mer", "15mer")

    # Test whether the seed.name is already defined. If not then use the input
    # start and stop coordinates
    if (seed.name %in% cannonical.seeds){
        target.seq = .get_canonical_seed(guide.seq, seed.name, allow_wobbles = allow_wobbles)

    } else if (!(seed.name %in% cannonical.seeds) & stop.pos & start.pos){
        target.seq = .get_custom_sequence(guide.seq, start.pos, stop.pos, allow_wobbles = allow_wobbles)

    } else {
        stop("Must provide either a canonical seed name or a custom seed name \
             with start.pos and stop.pos args.")

    }

    return(c(Seed.Name = seed.name,
             Guide = Biostrings::RNAString(guide.seq),
             Seed.Seq.RNA = Biostrings::RNAString(Biostrings::reverseComplement(target.seq)),
             Seed.Seq.DNA = Biostrings::reverseComplement(target.seq),
             Target.Seq = target.seq
             ))
}


#' Get the target sequence for canonical seed definitions
#'
#' @param guide.seq The input RNA string as a character vector. Must be at least
#' 8 characters long.
#' @param seed.name The name of the canonical seed definition to return.
#' Options: "mer7m8", "mer8", "mer6", "mer7A1", "Full", "18-mer", "15-mer"
#' @param allow_wobbles If true, all positions for a G:U wobble are considered.
#'
#' @return A `DNAString` of the canonical seed sequence that is the search term.
#' @export
#'
#' @examples
#' guide.seq = "UUAUAGAGCAAGAACACUGUUUU"
#'
#' canonical.seed.seq = SeedMatchR:::.get_canonical_seed(guide.seq,
#' seed.name = "mer7m8")
.get_canonical_seed <- function(guide.seq,
                                seed.name = c("mer7m8",
                                              "Full",
                                              "18mer",
                                              "15mer",
                                              "mer8",
                                              "mer6",
                                              "mer7A1"),
                                allow_wobbles = FALSE){
    seed.name = match.arg(seed.name)

    # Check parameters
    stopifnot(nchar(guide.seq) >= 8)

    if(allow_wobbles){
      guide.seq = stringr::str_replace_all(guide.seq, "U", "Y")
    }

    # Get the sub string of the sequence for the seed g2-g8
    base.seed = Biostrings::RNAString(substring(guide.seq, 2, 8))

    # Get the reverse complement of the sequence and convert to DNAString object
    seed.rev.comp.dna = Biostrings::reverseComplement(Biostrings::DNAString(base.seed))

    # Depending on seed.name, return the seed sequence: assume there will always
    # be an A at g1. This does not apply to the Full definition

    target.seq <- switch(
        seed.name,
        "Full" = .get_custom_sequence(guide.seq, 2, nchar(guide.seq), allow_wobbles = allow_wobbles),
        "18mer" = .get_custom_sequence(guide.seq, 2, 19, allow_wobbles = allow_wobbles),
        "15mer" = .get_custom_sequence(guide.seq, 2, 16, allow_wobbles = allow_wobbles),
        "mer8" = Biostrings::DNAString(paste0(substr(seed.rev.comp.dna, 1, 7),
                                              "A")),
        "mer7A1" = Biostrings::DNAString(paste0(substr(seed.rev.comp.dna, 2, 7),
                                                "A")),
        "mer7m8" = .get_custom_sequence(guide.seq, 2, 8, allow_wobbles = allow_wobbles),
        "mer6" = .get_custom_sequence(guide.seq, 2, 7, allow_wobbles = allow_wobbles)
    )

    return(target.seq)
}


#' Get a custom target sequence based on start and stop position
#'
#' This function will take a guide sequence and subset it based on custom start
#' and stop positions. The sequence will then be converted to the
#' reverse-complement `DNAString` of the input `RNAString` sequence. The start
#' and stop positions are 1-based indexes.
#'
#' @param guide.seq The input RNA string as a character vector
#' @param start.pos Start position index. 1-based index
#' @param stop.pos Stop position index.
#' @param allow_wobbles If true, all positions for a G:U wobble are considered.
#'
#' @return A `DNAString` of the customized search sequence
#' @export
#'
#' @examples
#' guide.seq = "UUAUAGAGCAAGAACACUGUUUU"
#'
#' custom.seq = SeedMatchR:::.get_custom_sequence(guide.seq, start.pos = 1,
#' stop.pos = 23)
.get_custom_sequence <- function(guide.seq, start.pos, stop.pos, allow_wobbles = FALSE){
    stopifnot(stop.pos <= nchar(guide.seq),
              start.pos > 0)

    if(allow_wobbles){
      guide.seq = stringr::str_replace_all(guide.seq, "U", "Y")
    }

    # Convert character vector to RNAString
    rna.string = Biostrings::RNAString(guide.seq)

    # Subset RNAString based on start and stop positions
    sub.seq = substring(rna.string, start.pos, stop.pos)

    seq.rna = Biostrings::RNAString(sub.seq)

    seq.dna = Biostrings::DNAString(seq.rna)

    target.seq = Biostrings::reverseComplement(seq.dna)

    return(target.seq)
}
