
################################################################################
# Tests for .get_canonical_seed
################################################################################
test_that("Extract canonical mer7m8", {
  expect_identical(SeedMatchR:::.get_canonical_seed(guide.seq = guide.seq,
                                                     seed.name = "mer7m8"),
                   Biostrings::DNAString("CTCTATA"))
})

test_that("Extract canonical mer8", {
  expect_identical(SeedMatchR:::.get_canonical_seed(guide.seq = guide.seq,
                                                    seed.name = "mer8"),
                   Biostrings::DNAString("CTCTATAA"))
})

test_that("Extract canonical mer7A1", {
  expect_identical(SeedMatchR:::.get_canonical_seed(guide.seq = guide.seq,
                                                    seed.name = "mer7A1"),
                   Biostrings::DNAString("TCTATAA"))
})

test_that("Extract canonical mer6", {
  expect_identical(SeedMatchR:::.get_canonical_seed(guide.seq = guide.seq,
                                                    seed.name = "mer6"),
                   Biostrings::DNAString("TCTATA"))
})

test_that("Extract Full sequence", {
  expect_identical(SeedMatchR:::.get_canonical_seed(guide.seq = guide.seq,
                                                    seed.name = "Full"),
                   Biostrings::DNAString("AAAACAGTGTTCTTGCTCTATA"))
})

test_that("Extract 18mer sequence", {
  expect_identical(SeedMatchR:::.get_canonical_seed(guide.seq = guide.seq,
                                                    seed.name = "18mer"),
                   Biostrings::DNAString("CAGTGTTCTTGCTCTATA"))
})
