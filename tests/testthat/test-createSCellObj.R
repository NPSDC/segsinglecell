context("test-input")

test_that("cDataPath var", {
    tsvFiles <- getTsvFiles("../../../research/out")
    segsList <- extSegs(tsvFiles)
    df <- createCountDf(segNames = segsList[['segNames']], tsvFiles = tsvFiles,
                        readEnd = segsList[['End']])
    annPath <- "../../../research/SraRunTable.tx"
    expect_that(createSCellObj(counts = , cDataPath = annPath, sepCData = ','),
                throws_error(paste(annPath, "is invalid path")))
})

context("test-output")
test_that("cDataPath var", {
    tsvFiles <- getTsvFiles("../../../research/out")
    segsList <- extSegs(tsvFiles)
    df <- createCountDf(segNames = segsList[['segNames']], tsvFiles = tsvFiles,
                        readEnd = segsList[['End']])
    annPath <- "../../../research/SraRunTable.txt"
    sc <- createSCellObj(counts = df, cDataPath = annPath, sepCData = ',')
    expect_equal(class(sc)[1], "SingleCellExperiment")
    expect_equal(ncol(colData(sc)), 34)
    expect_equal(nrow(colData(sc)), 35)

    tsvFiles <- getTsvFiles("../../../research/out_test/out_test_extSegs3/")
    segsList <- extSegs(tsvFiles)
    df <- createCountDf(segNames = segsList[['segNames']], tsvFiles = tsvFiles,
                        readEnd = segsList[['End']])
    annPath <- "../../../research/SraRunTable.txt"
    sc <- createSCellObj(counts = df, cDataPath = annPath, sepCData = ',')
    expect_equal(class(sc)[1], "SingleCellExperiment")
    expect_equal(ncol(colData(sc)), 34)
    expect_equal(nrow(colData(sc)), 2)
})
