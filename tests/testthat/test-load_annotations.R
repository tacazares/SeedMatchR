test_that("Expect correct filter for default settings.", {
  expect_snapshot(suppressMessages(SeedMatchR::build_annotation_filter(standard.chroms = TRUE,
                                                      remove.na.rows = TRUE,
                                                      protein.coding = TRUE))
  )
})

test_that("Expect correct filter for selecting symbol and entrez ID with TxBiotypeFilter enabled.", {
  expect_snapshot(suppressMessages(SeedMatchR::build_annotation_filter(standard.chroms = TRUE,
                                                      remove.na.rows = TRUE,
                                                      protein.coding = TRUE,
                                                      canonical = FALSE,
                                                      symbol = c("Ttr", "Hao1"),
                                                      entrez_id = c("1", "300"),
                                                      add.filter = AnnotationFilter::TxBiotypeFilter("Nonsense Mediated Decay")))
  )
})

test_that("Expect correct type for added AnnotationFilter.", {
  expect_snapshot(suppressMessages(SeedMatchR::build_annotation_filter(standard.chroms = TRUE,
                                                      remove.na.rows = TRUE,
                                                      protein.coding = TRUE,
                                                      canonical = TRUE,
                                                      symbol = c("Ttr", "Hao1"),
                                                      entrez_id = c("1", "300"),
                                                      add.filter = AnnotationFilter::TxBiotypeFilter("Nonsense Mediated Decay")))
  )
})

test_that("Expect error for missing value in added AnnotationFilter.", {
  expect_error(suppressMessages(SeedMatchR::build_annotation_filter(standard.chroms = TRUE,
                                                      remove.na.rows = TRUE,
                                                      protein.coding = TRUE,
                                                      canonical = TRUE,
                                                      symbol = c("Ttr", "Hao1"),
                                                      entrez_id = c("1", "300"),
                                                      add.filter = AnnotationFilter::TxBiotypeFilter()))
  )
})

test_that("Expect error when transcript support level is not less than or equal to 5.", {
  expect_error(suppressMessages(SeedMatchR::build_annotation_filter(standard.chroms = TRUE,
                                                      remove.na.rows = TRUE,
                                                      protein.coding = TRUE,
                                                      transcript.support = 6))
  )
})
