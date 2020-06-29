#' Gets the sample name(SRA) given a tsvFile which also contains the directory path
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

checkReadConsistency <- function(tsvFiles, cores = 1)
{
    library(data.table)
    library(doParallel)
    library(foreach)
    registerDoParallel(cores = cores)
    readEnds <- foreach(i = 1:length(tsvFiles), .inorder = F, .combine = 'c') %dopar%
    {
        tsvFile <- tsvFiles[i]
        getReadEnd(read.table(tsvFile, sep = "\t", stringsAsFactors = F, header = F, nrows = 1))
    }
    registerDoSEQ()
    gc()
    same <- T
    type <- readEnds[1]
    for(i in readEnds)
    {
      if(i != type)
      {
        same <- F
        break()
      }
    }
    if(!same)
      stop('Read Ends not consistent across files')
    return(type)
}

readSegCounts <- function(tsvFiles, segHashDf, cores = 1)
{
    library(data.table)
    library(doParallel)
    library(foreach)

    readEnd <- checkReadConsistency(tsvFiles)
    #cl <- makeCluster(cores)
    registerDoParallel(cores = cores)
    dfReq <- foreach(i = 1:length(tsvFiles), .inorder = F) %dopar%
    {
        tsvFile <- tsvFiles[i]
        sample <- getSample(tsvFile)
        df <- read.table(tsvFile, header = T, sep = "\t", stringsAsFactors = F)
        if(readEnd == "P")
            data.frame(Segs = paste(segHashDf[df$SEG1ID,], segHashDf[df$SEG2ID,], sep = '-'),
                   Sample = rep(sample, nrow(df)), Counts = df[,'count'], stringsAsFactors = F)
        else
            data.frame(Segs = segHashDf[df$SEG1ID,], Sample = rep(sample, nrow(df)),
                       Counts = df[,'count'], stringsAsFactors = F)
    }
    registerDoSEQ()
    gc()
    dfReq <- data.frame(rbindlist(dfReq))
    dfReq <- dfReq[order(dfReq$Segs),]
    return(dfReq)
}

createRowData <- function(gRange, rNamesCount)
{
    #rowAnnList <- values(unlist(rowAnnList[["gRange"]]))[,"segID","st","end"]
    readEnd <- if(grepl("[0-9]-[0-9]", rNamesCount[1])) "P" else "S"

    if(readEnd == "S")
    {
        ##To be completed
        gRListReq <- gRList[rNamesCount]
    }

    else
    {
        seg1 <- as.numeric(sapply(strsplit(rNamesCount, split = "-", fixed = T), function(x) x[1]))
        seg2 <- as.numeric(sapply(strsplit(rNamesCount, split = "-", fixed = T), function(x) x[2]))

        gRangeSeg1 <- gRange[seg1]
        gRangeSeg2 <- gRange[seg2]

        start1 <- start(gRangeSeg1)
        start2 <- start(gRangeSeg2)

        startReq <- start1
        startReq[start1 > start2] <- start2[start1 > start2]

        end1 <- end(gRangeSeg1)
        end2 <- end(gRangeSeg2)
        endReq <- end1
        endReq[end1 < end2] <- end2[end1 < end2]

        chr1 <- as.character(seqnames(gRangeSeg1))
        chr2 <- as.character(seqnames(gRangeSeg2))

        if(sum(chr1 != chr2) != 0)
            stop("chr1 not same as chr2")
        #chr <- paste(chr1, chr2, sep = '-')
        #print('ss')
        gRReq <- GRanges(seqnames = chr1, ranges = IRanges(startReq, endReq),
                      mcol = DataFrame(Seg1 = seg1, Seg2 = seg2))
    }
    return(gRReq)
}

createColData <- function(cPath, cNames, sep = ',')
{
    if(!(file.exists(cPath)))
      stop(paste(cPath, "is invalid path"))

    cData <- read.delim(cPath, sep = sep, stringsAsFactors = F, row.names = 1)

    # if(!('Run' %in% colnames(cData)))
    #   stop('Run not a column in cData')
    if(sum(!(cNames %in% rownames(cData))) != 0)
      stop("cNames not present")
    
    cDataReq <- cData[cNames,]
    return(cDataReq)
}


createSegAnnotation <- function(segsMetaFile, sep = '\t')
{
    library(IRanges)
    library(GenomicRanges)
    if(!(file.exists(segsMetaFile)))
        stop(paste(segsMetaFile, "is invalid path"))
    segMetaDf <- read.delim(file = segsMetaFile, header = T, sep = sep, stringsAsFactors = F)
    segMetaDf <- segMetaDf[order(segMetaDf[,"segID"]),]
    correctDfStarEnd <- function(df)
    {
      inds <- which(df$st > df$end)
      df[inds, c('st', 'end')] <- df[inds, c('end', 'st')]
      return(df)
    }

    segMetaDf <- correctDfStarEnd(segMetaDf)
    #segMetaDf$ind <- seq(nrow(segMetaDf))
    segHash <- data.frame(ind = seq(nrow(segMetaDf)), row.names = segMetaDf[,'segID'])

    if(sum(rownames(segHash) == segMetaDf[,'segID']) != nrow(segMetaDf))
        stop("rownames not match")
   # segHash <- segMetaDf[,c("segID", "ind")]
    SegGRanges <- makeGRangesFromDataFrame(segMetaDf, keep.extra.columns = T, start.field = "st",
                  end.field = "end", seqnames.field = "chrom", strand.field = "strand")
    #values(SegGRanges)[,"ind"] <- segHash[values(SegGRanges)[,'segID'], "ind"]
    remove(segMetaDf)
    gc()

    return(list("gRange" = SegGRanges, "hashDf" = segHash))
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

    sep = '/'
    if(endsWith(dir, '/'))
      sep = ''
    tsvFiles <- paste(dir, tsvFiles, sep = sep)
    return(tsvFiles)
}

#' Gets the type of read end
#'
#' Outputs whether a file is single or paired end read by checking whether arg starts with "SEG"
#' @param arg A character/number which is the 1st entry of the second Column of the TSV file
#'
#' @return character containing either 'S' or 'P' representing single or paired end read
#'
getReadEnd <- function(colNames)
{
    if(length(colNames) < 2)
        stop("ColNames length is less than 2")
    return(if(colNames[2] == "SEG2ID") "P" else "S")
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
createSCellObj <- function(dfReads, cData, rData, metaData)
{
    library(SingleCellExperiment)
    library(Matrix)

    dfReads[,"Segs"] <- as.factor(dfReads[,"Segs"])
    dfReads[,"Sample"] <- as.factor(dfReads[,"Sample"])
    mat <- sparseMatrix(as.numeric(dfReads[,"Segs"]), as.numeric(dfReads[,"Sample"]),
                      x = dfReads$Counts, dimnames = list(levels(dfReads$Segs),
                                                        levels(dfReads$Sample)), index1 = T)

    sce <- SingleCellExperiment(assays = list(counts = mat), colData = cData, rowData = rData, 
                                metadata = list("meta" = metaData))
    #int_elementMetadata(sce) <- DataFrame(rData)
    #int_metadata(sce) <- list(metaData)
    sce <- sce[,order(colnames(sce))]
    return(sce)
}

createSummExpObj <- function(dfReads, cData, rData, metaData)
{
    library(SummarizedExperiment)
}

#' Reads the segment counts from the directory containing tsv files
#'
#' @param dir character containing the directory in which tsv files are present
#' @param cDataPath character denoting path of the sample/cell annotation file
#' @param segMetaFile character denoting path of the seg meta annotation file
#' @param sepCData character denoting the separator to be used for reading
#' @param dataType S or B denoting single cell or bulk data
#' @param sep vector of length two containing separators for column Data file and seg Meta file
#'
#' @return SingleCellExperiment or SummarizedExperiment containing the counts
#'
#' @export
createSegCounts <- function(dir, cDataPath, segMetaPath, savePath, dataType = 'S', sep = c(",", "\t"), cores = 1)
{
    if(is.null(cDataPath) | !file.exists(cDataPath))
      stop("cDataPath is NULL or cDataPath does not exist")
    if(is.null(segMetaPath) | !file.exists(segMetaPath))
      stop("rDataPath is NULL or rDataPath does not exist")
    if(!(dataType %in% c("S", "B")))
      stop("dataType neither S nor B")
    if(length(sep) != 2)
      stop("Length of sep should be 2")

    tsvFiles <- getTsvFiles(dir)[1:10]
    print("Extracted Tsv files")

    segAnnotation <- createSegAnnotation(segsMetaFile = segMetaPath, sep = sep[2])
    print("Completed Seg Annotation")

    dfReads <- readSegCounts(tsvFiles = tsvFiles, segHashDf = segAnnotation[["hashDf"]], cores = cores)
    print("Extracted Reads")

    cData <- createColData(cPath = cDataPath, cNames = as.character(unique(dfReads[,"Sample"])))
    print("Extracted Column Information")

    rowGRanges <- createRowData(gRange = segAnnotation[["gRange"]],
                                rNamesCount = as.character(unique(dfReads[, "Segs"])))
    print("Extracted Row Information")

    if(dataType == "S")
        segCounts <- createSCellObj(dfReads = dfReads, cData = cData, rData = rowGRanges,
                            metaData = segAnnotation[["gRange"]])
    save(segCounts, file = savePath)
    return(segCounts)
}

