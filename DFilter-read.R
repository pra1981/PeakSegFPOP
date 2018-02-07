DFilter_peaks.vec <- Sys.glob("~/PeakSegFPOP/labels/H*/samples/*/*/DFilter_peaks.tsv")
N <- 100
for(DFilter_peaks.tsv in DFilter_peaks.vec){
  peaks.dt <- fread(paste("sed 's/ \\+/\t/g'", DFilter_peaks.tsv))[order(-`-log10(P-value)`, -maxScore)]
  sample.dir <- dirname(DFilter_peaks.tsv)
  print(sample.dir)
  fwrite(
    peaks.dt[, list(
      chromosome, `peak-start`, `peak-end`, mean=0)],
    file.path(sample.dir, "DFilter_peaks.bedGraph"),
    sep="\t", col.names=FALSE)
  topN.dt <- peaks.dt[1:N]
  fwrite(
    topN.dt[, list(
      chromosome, `peak-start`, `peak-end`, mean=0)],
    file.path(sample.dir, paste0("DFilter.top", N, "_peaks.bedGraph")),
    sep="\t", col.names=FALSE)
}
