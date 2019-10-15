context("test on output")

test_that("test output is S or P", {
    arg <- 1
    expect_equal(getReadEnd(arg), 'S')
    arg <- "SE"
    expect_equal(getReadEnd(arg), 'S')
    arg <- "SEG0001"
    expect_equal(getReadEnd(arg), 'P')
    arg <- "SEG1"
    expect_equal(getReadEnd(arg), 'P')
    arg <- "10"
    expect_equal(getReadEnd(arg), 'S')
})
