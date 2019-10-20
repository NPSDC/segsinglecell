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

#' Gets the segment count TSV files within a directory
#'
#' @param dir character containing the directory in which tsv files are expected
#'
#' @return Vector of the tsv files with the appropriate path
#'
getTsvFiles <- function(dir)
{
  if(!dir.exists(dir))
    stop(paste(dir, "is invalid directory"))

  tsvFiles <- list.files(dir, pattern = "*.tsv")
  if(length(tsvFiles) == 0)
    stop("No TSV files")

  tsvFiles <- tsvFiles[endsWith(tsvFiles, "tsv")]
  if(length(tsvFiles) == 0)
    stop("No TSV files")

  tsvFiles <- paste(dir, tsvFiles, sep = '/')
  return(tsvFiles)
}

#' Gets the type of read end
#'
#' Outputs whether a file is single or paired end read by checking whether arg starts with "SEG"
#' @param arg A character/number which is the 1st entry of the second Column of the TSV file
#'
#' @return character containing either 'S' or 'P' representing single or paired end read
#'
getReadEnd <- function(arg)
{
    read <- "S" #Single
    if(class(arg) == 'numeric' | class(arg) == 'integer')
      return("S")

    if(startsWith(arg, "SEG"))
      read <- "P" #Paired

    return(read)
}


#' Gets the combined segment pairs within the tsv files and the type of read end
#'
#' @param tsvFiles vector conatining the path of tsvFiles
#'
#' @return list with two arguments with first a character containing read end and
#'         second being another list containing segment/segment-pairs corresponding
#'         to each tsv file
#' @export
extSegs <- function(tsvFiles)
{
    segNames <- list()
    readEnd <- "S"
    for(i in seq_along(tsvFiles))
    {
        tsvFile <- tsvFiles[i]
        if(!file.exists(tsvFile))
            stop(paste(tsvFile, "is invalid path"))
        dfSegs <- read.table(tsvFile, header = F, sep = "\t", stringsAsFactors = F)

        if(i == 1)
            readEnd <- getReadEnd(dfSegs[1, "V2"])

        readEndCur <- getReadEnd(dfSegs[1, "V2"])

        if(readEndCur != readEnd)
            stop("All files are not single or paired end")

        if(readEnd == "S")
            segNames[[i]] <- as.character(dfSegs$V1)
        else
            segNames[[i]] <- paste(as.character(dfSegs$V1), as.character(dfSegs$V2), sep = "_")

        if(sum(duplicated(segNames[[i]])) != 0)
            stop(paste(tsvFiles[i], "Duplicated segments"))
    }
    return(list("End" = readEnd, "segNames" = segNames))
}

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

#' Creates a SingleCellExperiment object given count data and cell annotation file
#'
#' @param counts data.frame containing segment counts with rownames as segment/segment-pairs
#'               columns as cells
#' @param cDataPath character path containing the annotation file corresponding to cells
#' @param sepCData character containing the separator that would be used when reading annotation file
#'
#' @return SingleCellExperiment object
#'
#' @export
createSCellObj <- function(counts, cDataPath, sepCData = ",")
{
    library(SingleCellExperiment)
    if(!(file.exists(cDataPath)))
      stop(paste(cDataPath, "is invalid path"))

    cData <- read.delim(cDataPath, sep = sepCData, stringsAsFactors = F)

    if(!('Run' %in% colnames(cData)))
      stop('Run not a column in cData')

    commSampInds <- match(colnames(counts), as.character(cData$Run))
    cDataReq <- cData[commSampInds,]

    sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), colData = cDataReq)
    return(sce)
}

#' Reads the segment counts from the directory containing tsv files
#'
#' @param dir character containing the directory in which tsv files are present
#' @param sc boolean denoting whether SingleCellExperiment has to created
#' @param cDataPath character denoting path of the cell annotation file
#' @param sepCData character denoting the separator to be used for reading
#' cell annotation file
#'
#' @return data.frame or SingleCellExperiment containing the counts
#'
#' @export
readSegCounts <- function(dir, sc = F, cDataPath = NULL, sepCData = ",")
{
    tsvFiles <- getTsvFiles(dir)
    segInf <- extSegs(tsvFiles)
    if(sc)
    {
        if(is.null(cDataPath) | !file.exists(cDataPath))
          stop("cDataPath is NULL or cDataPath does not exist")
        else
        {
          countDf <- createCountDf(segNames = segInf[["segNames"]],
                                   tsvFiles = tsvFiles, readEnd = segInf[["End"]])
          countDf <- createSCellObj(counts = countDf, cDataPath = cDataPath, sepCData = sepCData)
          return(countDf)
        }
    }
    countDf <- createCountDf(segNames = segInf[["segNames"]],
                             tsvFiles = tsvFiles, readEnd = segInf[["End"]])
    return(countDf)
}

