context("test output")
test_that("out is df, nrows, ncols", {
    tsvFiles <- getTsvFiles("../../../research/out_test/out_test_extSegs3/")
    segsList <- extSegs(tsvFiles)

    df <- createCountDf(segNames = segsList[['segNames']], tsvFiles = tsvFiles,
                        readEnd = segsList[['End']])
    expect_equal(nrow(df), length(Reduce('union', segsList[['segNames']])))
    expect_equal(ncol(df), 2)

    tsvFiles <- getTsvFiles("../../../research/out")
    segsList <- extSegs(tsvFiles)

    df <- createCountDf(segNames = segsList[['segNames']], tsvFiles = tsvFiles,
                        readEnd = segsList[['End']])
    expect_equal(nrow(df), length(Reduce(union, segsList[['segNames']])))
    expect_equal(ncol(df), 35)
})
