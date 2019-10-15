context("test on inputs")

test_that("test for tsvFiles variable", {
    tsvFiles <- getTsvFiles("../../../research/out")
    tsvFilesCopy <- tsvFiles
    tsvFilesCopy[10] <- substr(tsvFiles[10], 1, length(tsvFiles[10]) - 1)
    expect_that(extSegs(tsvFilesCopy), throws_error(paste(tsvFilesCopy[10],
                                                    "is invalid path")))

    tsvFiles <- getTsvFiles("../../../research/out_test/out_test_extSegs/") #contains segment tsv files in which 2nd column has been removed in one of the files
    expect_that(extSegs(tsvFiles), throws_error("All files are not single or paired end"))

    tsvFiles <- getTsvFiles("../../../research/out_test/out_test_extSegs2//") #contains segment tsv files single end
    expect_that(extSegs(tsvFiles), throws_error(paste(tsvFiles[1], "Duplicated segments")))

})

context("test on output, length, each component further")
test_that("test output length, sub length, list, character, character expected", {
    tsvFiles <- getTsvFiles("../../../research/out")
    out <- extSegs(tsvFiles)

    expect_equal(length(out), 2)
    expect_equal(length(out[['End']]), 1)
    expect_equal(out[['End']], 'P')
    expect_equal(length(out[['segNames']]), 35)
    expect_equal(sum(sapply(out[['segNames']], function(x) sum(sapply(strsplit(x, split = '_', fixed = T),
                                                               length) == 2) == length(x)) == T), 35)
    tsvFiles <- getTsvFiles("../../../research/out_test/out_test_extSegs3") ##Single End reads with top 10 entries
    out <- extSegs(tsvFiles)
    expect_equal(length(out), 2)
    expect_equal(length(out[['End']]), 1)
    expect_equal(out[['End']], 'S')
    expect_equal(length(out[['segNames']]), 2)
})
