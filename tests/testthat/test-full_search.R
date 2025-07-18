test_that("Expect a full search for all default sequence definitions.", {
  expect_snapshot(SeedMatchR::full_search(guide.seq, seqs))
})
