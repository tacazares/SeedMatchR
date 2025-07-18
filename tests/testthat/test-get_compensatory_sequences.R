expected_cols <- c(
  "Seed_in_Transcript", "Compensatory_Region_in_Transcript", "TX_ID", "Distance",
  "Seed_Region_in_Guide", "Seed_Start", "Seed_End", "Compensatory_Start",
  "Compensatory_End", "Compensatory_Region_in_Guide", "mismatch_in_seed",
  "GC_Content", "CompensatoryLength", "relative_position_guide_seq_start",
  "relative_position_guide_seq_end"
)

test_that("get_compensatory_sequences returns expected columns and types", {
  guide.seq <- "AUGCUACGUGACUGAUGCUAGCUAGCUAGCUAGC"
  seed <- "AUGCUAC"
  sequences <- Biostrings::DNAStringSet(c("seq1" = "ATGCTACGTGACTGATGCTAGCTAGCTAGCTAGC",
                              "seq2" = "TTGCTACGTGACTGATGCTAGCTAGCTAGCTAGC"))

  result <- get_compensatory_sequences(
    guide.seq = guide.seq,
    seed = seed,
    sequences = sequences,
    compensatory_region_length_min = 6,
    compensatory_region_length_max = 8,
    limit = 100,
    mismatch_ = 1,
    mismatch_compensatory_region_ = 0
  )

  expect_true(is.data.frame(result))
  expect_true(all(expected_cols %in% colnames(result)))
})

test_that("get_compensatory_sequences handles no matches gracefully", {
  guide.seq <- "AUGCUACGUGACUGAUGCUAGCUAGCUAGCUAGC"
  seed <- "CCCCCCC"
  sequences <- Biostrings::DNAStringSet(c("seq1" = "ATGCTACGTGACTGATGCTAGCTAGCTAGCTAGC"))

  expect_message(
    result <- get_compensatory_sequences(
      guide.seq = guide.seq,
      seed = seed,
      sequences = sequences,
      compensatory_region_length_min = 6,
      compensatory_region_length_max = 8,
      limit = 100,
      mismatch_ = 0,
      mismatch_compensatory_region_ = 0
    ),
    "No compensatory regions found"
  )
  expect_true(is.data.frame(result))
  expect_true(all(expected_cols %in% colnames(result)))
  expect_equal(nrow(result), 0)
})

test_that("get_compensatory_sequences errors on invalid guide.seq", {
  guide.seq <- "AUGC"
  seed <- "AUGC"
  sequences <- Biostrings::DNAStringSet(c("seq1" = "ATGCTACGTGACTGATGCTAGCTAGCTAGCTAGC"))

  expect_error(
    get_compensatory_sequences(
      guide.seq = guide.seq,
      seed = seed,
      sequences = sequences
    ),
    "nchar\\(guide.seq\\) >= 8"
  )
})
