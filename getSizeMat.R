n.bins <- 20
n.bins.each.side <- n.bins/2
## Input a peaks_mat.csv file and output the coverage matrix.
readSizeMat <- function(peaks_mat.csv){
  mean.dt <- fread(peaks_mat.csv)
  mean.mat <- as.matrix(mean.dt[,-1,with=FALSE])
  rownames(mean.mat) <- mean.dt$V1
  colnames(mean.mat) <- NULL
  mean.mat
}
## Input a peaks.bedGraph file, and create a corresponding
## peaks_mat20000bp.csv file, a mean coverage matrix for n.bins of
## size bases.per.bin centered around the middle of the peak.
getSizeMat <- function(peaks.bedGraph, bases.per.bin=20000){
  sample.dir <- dirname(peaks.bedGraph)
  no.bedGraph <- sub(".bedGraph", "", basename(peaks.bedGraph))
  peaks_mat.csv <- file.path(
    sample.dir, paste0(
      no.bedGraph, "_mat",
      bases.per.bin, "bp",
      ".csv"))
  if(file.exists(peaks_mat.csv)){
    return(readSizeMat(peaks_mat.csv))
  }
  jpeaks <- fread(peaks.bedGraph)
  coverage.bigWig <- file.path(sample.dir, "coverage.bigWig")
  setnames(jpeaks, c("chrom", "peakStart", "peakEnd", "mean"))
  jpeaks[, peakBases := peakEnd-peakStart]
  jpeaks[, peakMid := as.integer((peakStart+peakEnd)/2)]
  bases.each.side <- n.bins.each.side * bases.per.bin
  list.dt <- jpeaks[, {
    from <- peakMid-bases.each.side
    to <- peakMid+bases.each.side
    zoom.coverage <- readBigWig(coverage.bigWig, chrom, max(0, from), to)
    edge.vec <- seq(from, to, l=n.bins+1)
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
