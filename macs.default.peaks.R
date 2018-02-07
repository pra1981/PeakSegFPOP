library(data.table)
peak.file.vec <- c(
  Sys.glob("~/peaks/*/*/*/*c=1.30103.narrowPeak"),
  Sys.glob("~/peaks/*/*/*c=1.30103.narrowPeak"))
peaks.dt <- data.table(peak.file=peak.file.vec)[, {
  sample.dir <- dirname(peak.file)
  maybe.exp.dir <- dirname(sample.dir)
  maybe.model.dir <- dirname(maybe.exp.dir)
  if(basename(maybe.model.dir)=="macs2.broad"){
    experiment <- basename(maybe.exp.dir)
    model <- "macs2.broad"
  }else{
    experiment <- if(basename(maybe.exp.dir)=="macs2.new.34c")"H3K4me3" else "H3K36me3"
    model <- "macs2.default"
  }
  sample.id <- basename(sample.dir)
  dt <- fread(peak.file)
  setnames(dt, c("chrom", "peakStart", "peakEnd"))
  data.table(experiment, model, sample.id, dt)
}, by=list(peak.file)]
macs.default.peaks <- peaks.dt[, -1, with=FALSE]
macs.default.peaks[, table(sample.id, experiment, model)]

fwrite(macs.default.peaks, "macs.default.peaks.csv")
