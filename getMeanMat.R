n <- 5
## Input a peaks_mat.csv file and output the coverage matrix.
readMeanMat <- function(peaks_mat.csv){
  mean.dt <- fread(peaks_mat.csv)
  mean.mat <- as.matrix(mean.dt[,-1,with=FALSE])
  rownames(mean.mat) <- mean.dt$V1
  colnames(mean.mat) <- NULL
  mean.mat
}
## Input a peaks.bedGraph file, and create a corresponding
## peaks_mat.csv file, a mean coverage matrix that we will use to make
## the meta-peak plot.
getMeanMat <- function(peaks.bedGraph){
  sample.dir <- dirname(peaks.bedGraph)
  no.bedGraph <- sub(".bedGraph", "", basename(peaks.bedGraph))
  peaks_mat.csv <- file.path(sample.dir, paste0(no.bedGraph, "_mat.csv"))
  if(file.exists(peaks_mat.csv)){
    return(readMeanMat(peaks_mat.csv))
  }
  jpeaks <- fread(peaks.bedGraph)
  coverage.bigWig <- file.path(sample.dir, "coverage.bigWig")
  setnames(jpeaks, c("chrom", "peakStart", "peakEnd", "mean"))
  list.dt <- jpeaks[, {
    zoomBases <- as.integer((peakEnd-peakStart)/2)
    from <- peakStart-zoomBases
    to <- peakEnd+zoomBases
    zoom.coverage <- readBigWig(coverage.bigWig, chrom, max(0, from), to)
    edge.vec <- c(
      seq(from, peakStart, l=n+1),
      seq(peakStart, peakEnd, l=2*n+1)[-1],
      seq(peakEnd, to, l=n+1)[-1])
    n.edges <- length(edge.vec)
    box.dt <- data.table(boxStart=edge.vec[-n.edges], boxEnd=edge.vec[-1])
    box.dt[, box.i := 1:.N]
    setkey(box.dt, boxStart, boxEnd)
    zoom.coverage[, chromStart1 := chromStart+1]
    setkey(zoom.coverage, chromStart1, chromEnd)
    over.dt <- foverlaps(zoom.coverage, box.dt, nomatch=0L)
    over.dt[, start := ifelse(chromStart1 < boxStart, boxStart, chromStart1)]
    over.dt[, end := ifelse(boxEnd < chromEnd, boxEnd, chromEnd)]
    over.dt[, weight := end-start]
    mean.dt <- over.dt[, list(
      mean=sum(weight*count)/(boxEnd-boxStart)
      ), by=list(box.i, boxStart, boxEnd)]
    box.dt[, mean := 0]
    box.dt[mean.dt$box.i, mean := mean.dt$mean]
    data.table(means=list(box.dt$mean))
  }, by=list(chrom, peakStart, peakEnd)]
  mean.mat <- do.call(rbind, list.dt$means)
  rownames(mean.mat) <- list.dt[, sprintf("%s:%d-%d", chrom, peakStart, peakEnd)]
  write.table(mean.mat, peaks_mat.csv, col.names=FALSE, quote=FALSE)
  mean.mat
}
