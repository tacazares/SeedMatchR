
test_that("Expect printed output from running function", {
  expect_snapshot(SeedMatchR::deseq_fc_ecdf(res,
                                         list("Background" = background.list,
                                              "mer7m8" = mer7m8.list))$stats)
})

test_that("Expect error because one list is empty", {
  expect_error(SeedMatchR::deseq_fc_ecdf(res,
                                          list("Background" = c(NULL),
                                               "mer7m8" = mer7m8.list))$stats)
})

test_that("Expect a change in line color", {
  expect_snapshot(SeedMatchR::deseq_fc_ecdf(res,
                                         list("Background" = background.list,
                                              "mer7m8" = mer7m8.list),
                                         palette = c("red", "green"))$stats)
})


test_that("Expect a change in plot ordering of lines", {
  expect_snapshot(SeedMatchR::deseq_fc_ecdf(res,
                                          list("Background" = background.list,
                                               "mer7m8" = mer7m8.list),
                                          palette = c("red", "green"),
                                          factor_order = c("mer7m8", "Background"))$stats)
})

test_that("Test parameters stats_test = Wilcoxen, alternative = less.",{
  expect_snapshot(SeedMatchR::ecdf_stat_test(res,
                                             background.list,
                                             mer7m8.list,
                                             stats_test = "Wilcoxen",
                                             alternative = "less"))
})

test_that("Check Wass statistic calculation.", {
  expect_snapshot(SeedMatchR::ecdf_stat_test(res,
                                          background.list,
                                          mer7m8.list,
                                          stats_test = "Wass")$statistic)
})
