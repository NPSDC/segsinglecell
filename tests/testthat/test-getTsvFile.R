context("test on inputs")

test_that("test for dir variable", {
    dir <- "out"
    expect_that(getTsvFiles(dir), throws_error(paste(dir, "is invalid directory")))

    dir <- "../../../research/"
    expect_that(getTsvFiles(dir), throws_error(paste("No TSV files")))

    dir <- "../../../research/out_test/out_test_getTsv/"
    expect_that(getTsvFiles(dir), throws_error(paste("No TSV files")))
})

context("test on outputs")
test_that("test output is a character, length", {
    dir <- "../../../research/out"
    expect_that(getTsvFiles(dir), is_a('character'))
    expect_equal(length(getTsvFiles(dir)), 35)
})
