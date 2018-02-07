library(data.table)
peak.file.vec <- c(
  Sys.glob("~/PeakSegFPOP/labels/H3K4me3_TDH_immune/samples/*/*/macs.*/peaks.csv"),
  Sys.glob("~/PeakSegFPOP/labels/H3K36me3_TDH_immune/samples/*/*/macs.*/peaks.csv"))

peaks.dt <- data.table(peak.file=peak.file.vec)[, {
  model.dir <- dirname(peak.file)
  model <- basename(model.dir)
  experiment <- sub(".*/labels/", "", sub("_.*", "", peak.file))
  sample.dir <- dirname(model.dir)
  sample.id <- basename(sample.dir)
  dt <- fread(peak.file, select=c(1:3, 6))
  data.table(experiment, model, sample.id, dt)
}, by=list(peak.file)]

macs.pvalues <- peaks.dt[, -1, with=FALSE]
macs.pvalues[, table(sample.id, experiment, model)]

fwrite(macs.pvalues, "macs.pvalues.csv")
