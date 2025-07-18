#' A modified vmatchPattern function that supports indels and mismatches
#'
#' @param pattern String to search for
#' @param subject DNAStringSet to query
#' @param max.mismatch Maximum number of mismatches allowed
#' @param min.mismatch Minimum number of mismatches allowed
#' @param with.indels Boolean for whether to include indels
#' @param fixed Boolean for whether to match character symbols
#' @param algorithm Which search algorithm to use. Defaults to "auto"
#'
#' @return A data frame of interval ranges per sequence
#' @export
#'
#' @examplesIf interactive()
#' print("TBD")
vmatchPattern2 <- function(pattern,
                           subject,
                           max.mismatch = 0,
                           min.mismatch = 0,
                           with.indels = FALSE,
                           fixed = TRUE,
                           algorithm = "auto")
{
  pattern <- Biostrings::DNAString(pattern)

  max.mismatch <- as.integer(max.mismatch)
  min.mismatch <- as.integer(min.mismatch)

  fixed <- .normalize_fixed_arg(fixed, subject)

  algo = .get_algo("auto", pattern,
                   max.mismatch, min.mismatch,
                   with.indels, fixed)

  C_ans <- S4Vectors::.Call2(
    "XStringSet_vmatch_pattern",
    pattern,
    subject,
    max.mismatch,
    min.mismatch,
    with.indels,
    fixed,
    "indels",
    "MATCHES_AS_RANGES",
    PACKAGE = "Biostrings"
  )

  unlisted_ans <- IRanges::IRanges(
    start = unlist(C_ans[[1L]],
                   use.names = FALSE),
    width = unlist(C_ans[[2L]],
                   use.names = FALSE)
  )

  ir.list = IRanges::relist(unlisted_ans, C_ans[[1L]])

  names(ir.list) <- names(subject)

  ir.list <- ir.list[!sapply(ir.list, function(x) length(x) == 0)]

  ir.df = unlist(ir.list) %>% as.data.frame()

  return(ir.df)
}

.normalize_fixed_arg <- function(fixed, subject, argname = "fixed")
{
  if (!is.logical(fixed) && !is.character(fixed))
    stop("'", argname, "' not a logical or character vector")
  if (is.logical(fixed)) {
    if (any(is.na(fixed)))
      stop("'", argname, "' has NAs")
    fixed_names <- names(fixed)
    if (is.null(fixed_names)) {
      if (!(length(fixed) %in% 1:2))
        stop("when an unnamed logical vector, '", argname,
             "' fixed must be of length 1 or 2")
      if (length(fixed) == 1)
        fixed <- c(fixed, fixed)
    }
    else {
      if (length(fixed) != 2)
        stop("when a named logical vector, '", argname,
             "' must be of length 2")
      if (!setequal(fixed_names, c("pattern", "subject")))
        stop("'", argname, "' names must be \"pattern\" and \"subject\"")
      fixed <- c(fixed["pattern"], fixed["subject"])
    }
  }
  else if (is.character(fixed)) {
    if (any(duplicated(fixed)) || !all(fixed %in% c("pattern",
                                                    "subject")))
      stop("when a character vector, '", argname, "' must be ",
           "a subset of 'c(\"pattern\", \"subject\")' ",
           "with no duplicated")
    fixed <- c("pattern" %in% fixed, "subject" %in% fixed)
  }
  if (!all(fixed) && !(Biostrings::seqtype(subject) %in% c("DNA", "RNA")))
    stop("'", argname, "' value only supported for a DNA or RNA subject ",
         "(you can only use 'fixed=TRUE' with your subject)")
  fixed
}

.get_algo <- function(algo, pattern, max.mismatch, min.mismatch, with.indels,
                     fixed)
{
  algos <- .valid_algos(pattern, max.mismatch, min.mismatch,
                        with.indels, fixed)
  if (algo == "auto")
    return(algos[1])
  if (!(algo %in% algos))
    stop(
      "valid algos for your string matching problem (best suited first): ",
      paste(paste("\"", algos, "\"", sep = ""), collapse = ", ")
    )
  algo
}

.valid_algos <- function(pattern, max.mismatch, min.mismatch, with.indels, fixed){
  {
    if (methods::is(pattern, "XString")) {
      pattern_min_length <- pattern_max_length <- length(pattern)
    }
    else if (methods::is(pattern, "XStringSet")) {
      pattern_min_length <- min(BiocGenerics::width(pattern))
      pattern_max_length <- max(BiocGenerics::width(pattern))
    }
    else {
      stop("'pattern' not an XString or XStringSet object")
    }
    if (pattern_min_length == 0L)
      stop("empty patterns are not supported")
    if (pattern_max_length > 20000L)
      stop("patterns with more than 20000 letters are not supported")
    if (max.mismatch != 0L && with.indels) {
      if (min.mismatch != 0L)
        stop("'min.mismatch' must be 0 when 'with.indels' is TRUE")
      return("indels")
    }
    algos <- character(0)
    if (max.mismatch == 0L && all(fixed)) {
      algos <- c(algos, "boyer-moore")
      if (pattern_max_length <= 64)
        algos <- c(algos, "shift-or")
      algos <- c(algos, "naive-exact")
    }
    else {
      if (min.mismatch == 0L && fixed[1] == fixed[2] && pattern_max_length <=
          64)
        algos <- c(algos, "shift-or")
    }
    c(algos, "naive-inexact")
  }
}
