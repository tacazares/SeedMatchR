test_that("Generate default report table.", {
  expect_snapshot(
    SeedMatchR::SeedMatchReport(
      seqs = seqs,
      res = res,
      guide.seq = guide.seq
    )$report
  )
})
#
# test_that("Generate indel report table.", {
#   expect_snapshot(
#     SeedMatchR::SeedMatchReport(
#       seqs = seqs,
#       res = res,
#       guide.seq = guide.seq,
#       indel.bool = TRUE
#     )$report
#   )
# })
#
# test_that("Generate wobble report table.", {
#   expect_snapshot(
#     SeedMatchR::SeedMatchReport(
#       seqs = seqs,
#       res = res,
#       guide.seq = guide.seq,
#       allow_wobbles = TRUE
#     )$report
#   )
# })
#
# test_that("Generate wobble + indel report table.", {
#   expect_snapshot(
#     SeedMatchR::SeedMatchReport(
#       seqs = seqs,
#       res = res,
#       guide.seq = guide.seq,
#       allow_wobbles = TRUE,
#       indel.bool = TRUE
#     )$report
#   )
# })
