
test_that("Expect correct counts for DESEQ2 mode for mer8", {
  expect_snapshot(SeedMatchR::SeedMatchR(res = res,
                                         seqs = seqs,
                                         sequence = guide.seq,
                                         seed.name = "mer7m8",
                                         res.format = "DESEQ2")
                  )
})

test_that("Expect correct counts for DESEQ2 mode for mer8", {
  expect_snapshot(SeedMatchR::SeedMatchR(res = res,
                                         seqs = seqs,
                                         sequence = guide.seq,
                                         seed.name = "mer8",
                                         res.format = "DESEQ2")
  )
})

test_that("Expect start and stop positions for DESEQ2 mode for custom sequence", {
  expect_error(SeedMatchR::SeedMatchR(res = res,
                                         seqs = seqs,
                                         sequence = guide.seq,
                                         seed.name = "Custom",
                                         res.format = "DESEQ2")
  )
})

test_that("Expect correct matches for granges mode for mer7m8", {
  expect_snapshot(SeedMatchR::SeedMatchR(seqs = seqs,
                                         sequence = guide.seq,
                                         seed.name = "mer7m8",
                                         res.format = "granges")
  )
})

test_that("Expect correct matches for iranges mode for mer7m8", {
  expect_snapshot(SeedMatchR::SeedMatchR(seqs = seqs,
                                         sequence = guide.seq,
                                         seed.name = "mer7m8",
                                         res.format = "iranges")
  )
})

test_that("Expect correct matches for data frame mode for mer7m8 with G:U wobbles", {
  expect_snapshot(SeedMatchR::SeedMatchR(seqs = seqs,
                                         sequence = guide.seq,
                                         seed.name = "mer7m8",
                                         res.format = "data.frame",
                                         allow_wobbles = TRUE,
                                         fixed = FALSE)
  )
})
