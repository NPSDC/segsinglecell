createSparseMatrix2 <- function(tsvFiles, segNames)
{
    library(Matrix)
    #library(forcats)
    #segs <- extSegs(tsvFiles = tsvFiles)
    print(length(segNames))
    print(length(tsvFiles))
    mat <- sparseMatrix(i = length(segNames), j = length(tsvFiles), x = 0,
                        dimnames = list(segNames, sapply(tsvFiles, getSample)))

    diag(mat) <- 0

    for(tsvFile in tsvFiles)
    {
        print(tsvFile)
        samp <- getSample(tsvFile)
        df <- read.table(tsvFile, header = T, sep = "\t", stringsAsFactors = F)
        segPairs <- paste(df$SEG1ID, df$SEG2ID, sep = '-')
        mat[segPairs, samp] <- df$count
    }
    return(mat)
}

#' Creates a binary Sparse Matrix given segment Meta file
createBinaryInitSparseMatrix <- function(metaFile)
{
    library(Matrix)
    if(!file.exists(metaFile))
        stop(paste(metaFile, "is invalid Path"))

    end <- 'PE'

    segsMetaDf <- read.delim(metaFile, stringsAsFactors = F)
    segs <- segsMetaDf[,1]

    if(end == 'SE'){
    }

    else
    {
        sparseInit <- sparseMatrix(1:length(segs), 1:length(segs), dimnames = list(segs, NULL))
        diag(sparseInit) <- F
    }
    return(sparseInit)
}

createBinarySparseMatrix <- function(tsvFiles, metaFile)
{
    sparseBinaryMat <- createBinaryInitSparseMatrix(metaFile = metaFile)
    for(tsvFile in tsvFiles)
    {
        dfSegs <- read.table(file = tsvFile, sep="\t", header = F,
                             colClasses = c(NA, NA, "NULL", "NULL", "NULL", "NULL", "NULL",
                                            "NULL", "NULL", "NULL"), stringsAsFactors = F)
        #print(rownames(sparseBinaryMat))
        rowInds <- match(dfSegs$V1, sparseBinaryMat@Dimnames[[1]])
        colInds <- match(dfSegs$V2, sparseBinaryMat@Dimnames[[1]])

        oInds <- order(rowInds)
        rowInds <- sort(rowInds)
        tableInds <- table(rowInds)

        if(sum(is.na(rowInds) != 0))
            stop("Segment is row is missing")
        if(sum(is.na(colInds) != 0))
            stop("Segment is col is missing")

        start = 1
        for(i in seq_along(as.numeric(names(tableInds))))
        {
            end = start + tableInds[i] - 1
            sparseBinaryMat[as.numeric(names(tableInds)[i]), colInds[oInds[start:end]]] = T
            start = end + 1
        }

        print('sup')
    }
    return(sparseBinaryMat)
}

#' Gets the combined segment pairs within the tsv files and the type of read end
#'
#' @param tsvFiles vector conatining the path of tsvFiles
#'
#' @return list with two arguments with first a character containing read end and
#'         second being another list containing segment/segment-pairs corresponding
#'         to each tsv file
#' @export
# extSegs <- function(tsvFiles, union = T, length = 100, cores = 1)
# {
#     segNamesList <- mclapply(tsvFiles, function(tsvFile)
#     {
#         segDf <- read.table(tsvFiles, header = T, sep = '\t')
#         segnames <- paste()
#     })
#     if(union)
#       segNamesUnion <- c()
#     else
#
#
# }

#' Creates a countMatrix corresponding to the tsv files
#'
#' @param segNames list containing the segment/segment-pairs corresponding to each tsv file
#' @param tsvFiles vector conatining the path of tsvFiles
#' @param readEnd character containing read end
#'
#' @return data.frame containing segment counts with rownames as segment/segment-pairs
#'         columns as cells
createCountDf <- function(segNames, tsvFiles, readEnd)
{
    ##Contains the union of all the segments present across files
    allSegs <- Reduce(union, segNames)
    cellNames <- sapply(strsplit(sapply(strsplit(tsvFiles, split = "/", fixed = T),
                                        function(x) x[[length(x)]]), split = '.',
                                 fixed = T), function(x) x[[1]])

    countDf <- data.frame(matrix(vector(), length(allSegs), length(tsvFiles),
                                 dimnames=list(allSegs, cellNames)))
    countDf[,] <- 0

    for(i in seq_along(tsvFiles))
    {
        tsvFile <- tsvFiles[i]
        df <- read.table(tsvFile, header = F, sep = "\t", stringsAsFactors = F)

        if(readEnd == "S")
            countDf[segNames[[i]],i] <- rowSums(df[,-1])
        else
            countDf[segNames[[i]],i] <- rowSums(df[,-c(1,2)])
    }
    return(countDf)
}

#' Creates the sparse Matrix from tsvFiles
createSparseMatrix <- function(tsvFiles, end = "PE", cores = 1)
{
    library(Matrix)
    library(data.table)
    #library(Parallel)
    library(doParallel)
    library(foreach)
    #dfReq <- data.frame(Segs = character(), Sample = character(), Counts = integer())

    # dfList <- mclapply(tsvFiles, function(tsvFile)
    # {
    #     sample <- getSample(tsvFile)
    #     df <- read.table(tsvFile, header = T, sep = "\t", stringsAsFactors = F)
    #     data.frame(Segs = paste(df$SEG1ID, df$SEG2ID, sep = '-'),
    #                Sample = rep(sample, nrow(df)), Counts = df[,'count'])
    #
    #     # dfReq <- rbind(dfReq, data.frame(Segs = paste(df$V1, df$V2, sep = '-'),
    #     #                                  Sample = rep(sample, nrow(df)),
    #     #                                  Counts = rowSums(df[,3:ncol(df)])))
    # }, mc.cores = cores)
    #dfReq <- transform(dfReq, Segs = factor(Segs), Sample  = factor(Sample))
    #  dfReq <- rbindlist(dfList)

    registerDoParallel(cores = cores)
    dfReq <- foreach(i = 1:length(tsvFiles), .inorder = F) %dopar%
        {
            tsvFile <- tsvFiles[i]
            sample <- getSample(tsvFile)
            df <- read.table(tsvFile, header = T, sep = "\t", stringsAsFactors = F)
            data.frame(Segs = paste(df$SEG1ID, df$SEG2ID, sep = '-'),
                       Sample = rep(sample, nrow(df)),
                       Counts = df[,'count'])
        }
    dfReq <- rbindlist(dfReq)
    registerDoSEQ()
    #return(dfReq)

    gc()
    dfReq <- dfReq[order(dfReq$Segs),]
    # mat <- sparseMatrix(as.numeric(dfReq$Segs), as.numeric(dfReq$Sample),
    #                     x = dfReq$Counts, dimnames = list(levels(dfReq$Segs),
    #                                                       levels(dfReq$Sample)))
    #
    # return(mat=mat)
}
