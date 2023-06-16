# test_that("You can retreive the correct mer7m8 seed sequence", {
#   expect_equal(get.seed("UUAUAGAGCAAGAACACUGUUUU", "mer7m8"), "CTCTATA")
# })
#
# test_that("You can retreive the correct mer7A1 seed sequence", {
#   expect_equal(get.seed("UUAUAGAGCAAGAACACUGUUUU", "mer7A1"), "TCTATAA")
# })
#
# test_that("You can retreive the correct mer8 seed sequence", {
#   expect_equal(get.seed("UUAUAGAGCAAGAACACUGUUUU", "mer8"), "CTCTATAA")
# })
#
# test_that("You can retreive the correct mer6 seed sequence", {
#   expect_equal(get.seed("UUAUAGAGCAAGAACACUGUUUU", "mer6"), "TCTATA")
# })
#
# test_that("get.seed() errors if input length < 8", {
#   expect_error(get.seed("UUAUAG", "mer6"))
# })
#
# test_that("get.seed() errors if input seq is not a character", {
#   expect_error(get.seed(1, "mer6"))
# })
#
# test_that("get.seed() errors if input seed.name is not a character", {
#   expect_error(get.seed("UUAUAGAGCAAGAACACUGUUUU", 1))
# })
#
# test_that("get.seed() errors if input seed.name is not in defined list", {
#   expect_error(get.seed("UUAUAG", "mer5"))
# })
