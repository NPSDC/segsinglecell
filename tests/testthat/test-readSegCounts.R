context("test-input")

test_that("cDataPath var", {
    dir <- "../../../research/out"
    expect_that(readSegCounts(dir = dir, sc = T, cDataPath = NULL, sepCData = ','),
                throws_error("cDataPath is NULL"))
    expect_that(readSegCounts(dir = dir, sc = T, cDataPath = '../../../research/SraRunTable.txt', sepCData = ','),
                throws_error("cDataPath is NULL"))
})

context("test-output")

test_that("cDataPath var", {
    dir <- "../../../research/out"
    df <- readSegCounts(dir = dir)
    expect_equal(class(df), "data.frame")

    df <- readSegCounts(dir = dir, sc = T, cDataPath = "../../../research/SraRunTable.txt")
    expect_equal(class(df)[1], "SingleCellExperiment")
})
