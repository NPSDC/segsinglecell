getSample <- function(sample)
{
    if(length(sample) != 1)
      stop("Length of file should be 1")
    if(grepl('/', sample))
    {
      sample <- strsplit(sample, split = "/", fixed = T)
      sample <- sample[[1]][length(sample[[1]])]
    }

    sample <- strsplit(sample, split = ".", fixed = T)
    return(sample[[1]][1])
}

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
    mat <- sparseMatrix(as.numeric(dfReq$Segs), as.numeric(dfReq$Sample),
                         x = dfReq$Counts, dimnames = list(levels(dfReq$Segs),
                                                           levels(dfReq$Sample)))

    return(mat=mat)
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

