library(PeakError)
library(data.table)
library(PeakSegPipeline)
library(ggplot2)
macs.pvalues <- fread("macs.pvalues.csv")
macs.pvalues[, peakBases := end - start]
setkey(macs.pvalues, experiment, sample.id)

## test AUC.
lik.gz.vec <- grep(
  "immune|other",
  Sys.glob("labels/H*_TDH_*/peaks_matrix_likelihood.tsv.gz"),
  value=TRUE)
lik.dt.list <- list()
for(lik.i in seq_along(lik.gz.vec)){
  lik.gz <- lik.gz.vec[[lik.i]]
  cat(sprintf("%4d / %4d %s\n", lik.i, length(lik.gz.vec), lik.gz))
  lik.dt <- fread(paste("zcat", lik.gz))
  pattern <- paste0(
    "(?<chrom>chr[^:]+)",
    ":",
    "(?<peakStart>[^-]+)",
    "-",
    "(?<peakEnd>.*)")
  pos.dt <- data.table(
    namedCapture::str_match_named(lik.dt$peak.name, pattern, list(
      peakStart=as.integer,
      peakEnd=as.integer)))
  set.name <- basename(dirname(lik.gz))
  lik.dt.list[[set.name]] <- data.table(pos.dt, lik.dt)
}

labels.bed.vec <- grep(
  "immune|other",
  Sys.glob("labels/H*_TDH_*/samples/*/*/labels.bed"),
  value=TRUE)
auc.dt.list <- list()
for(labels.i in seq_along(labels.bed.vec)){
  labels.bed <- labels.bed.vec[[labels.i]]
  cat(sprintf("%4d / %4d %s\n", labels.i, length(labels.bed.vec), labels.bed))
  labels.dt <- fread(labels.bed)
  setnames(labels.dt, c("chrom", "chromStart", "chromEnd", "annotation"))
  sample.dir <- dirname(labels.bed)
  sample.id <- basename(sample.dir)
  group.dir <- dirname(sample.dir)
  samples.dir <- dirname(group.dir)
  proj.dir <- dirname(samples.dir)
  set.name <- basename(proj.dir)
  experiment <- sub("_.*", "", set.name)
  test.samples <- sub(".*_", "", set.name)
  train.samples <- ifelse(test.samples=="immune", "other", "immune")
  train.proj.dir <- paste0(experiment, "_TDH_", train.samples)
  sample.path <- paste0(basename(group.dir), "/", sample.id)
  cols.vec <- c("chrom", "peakStart", "peakEnd", sample.path)
  joint.peaks.dt <- lik.dt.list[[train.proj.dir]][, cols.vec, with=FALSE]
  setnames(joint.peaks.dt, c("chrom", "peakStart", "peakEnd", "lik"))
  setkey(joint.peaks.dt, lik)
  select.dt <- data.table(experiment, sample.id)
  macs.peaks.dt <- macs.pvalues[select.dt, list(
    model, chrom=chr, peakStart=start, peakEnd=end, lik=`-log10(pvalue)`)]
  pred.dt <- rbind(
    data.table(model="PeakSeg", joint.peaks.dt),
    macs.peaks.dt)
  setkey(labels.dt, chrom, chromStart, chromEnd)
  setkey(pred.dt, chrom, peakStart, peakEnd)
  over.dt <- foverlaps(pred.dt, labels.dt, nomatch=0L)
  thresh.diff.dt <- over.dt[, {
    ord.peaks <- .SD[order(-lik)]
    one.label <- data.table(chromStart, chromEnd, annotation)
    thresh.dt <- data.table(thresh.i=0:.N)[, {
      peaks.df <- if(thresh.i==0){
        Peaks()
      }else{
        ord.peaks[1:thresh.i, list(chromStart=peakStart, chromEnd=peakEnd)]
      }
      err.df <- PeakErrorChrom(peaks.df, one.label)
      data.table(thresh=ifelse(thresh.i==0, Inf, ord.peaks$lik[thresh.i]), err.df)
    }, by=list(thresh.i)]
    thresh.dt[, list(thresh=thresh[-1], diff.tp=diff(tp), diff.fp=diff(fp))]
  }, by=list(model, chrom, chromStart, chromEnd, annotation)]
  joint.err <- PeakError(
    joint.peaks.dt[, list(chrom, chromStart=peakStart, chromEnd=peakEnd)],
    labels.dt)
  total.dt <- with(joint.err, data.table(
    tp=sum(tp),
    fp=sum(fp),
    possible.tp=sum(possible.tp),
    possible.fp=sum(possible.fp)))
  roc.dt <- thresh.diff.dt[, {
    thresh.ord <- .SD[, list(
      diff.tp=sum(diff.tp),
      diff.fp=sum(diff.fp)
    ), by=list(thresh)][order(-thresh)]
    some.roc <- thresh.ord[, data.table(
      thresh=c(Inf, thresh),
      tp=cumsum(c(0, diff.tp)),
      fp=cumsum(c(0, diff.fp))
      )]
    all.roc <- rbind(
      some.roc,
      data.table(
        thresh=c(0, -1),
        tp=c(some.roc[.N, tp], 0),
        fp=total.dt$possible.fp))
    data.table(all.roc, total.dt[, list(possible.fp, possible.tp)])
  }, by=list(model)]
  roc.dt[, FPR := fp/possible.fp]
  roc.dt[, TPR := tp/possible.tp]
  roc.total <- roc.dt[model=="PeakSeg" & 0 < thresh][which.min(thresh)]
  stopifnot(roc.total$tp==total.dt$tp)
  stopifnot(roc.total$fp==total.dt$fp)
  auc.dt <- roc.dt[, list(auc=geometry::polyarea(FPR, TPR)), by=model]
  print(auc.dt)
  auc.dt.list[[paste(experiment, sample.id)]] <- data.table(
    experiment, sample.id, test.samples, auc.dt)
}
auc.dt <- do.call(rbind, auc.dt.list)

auc.wide <- dcast(auc.dt, experiment + sample.id + test.samples ~ model)
ggplot()+
  geom_abline(slope=1, intercept=0, color="grey")+
  geom_point(aes(
    PeakSeg, macs.default),
             shape=1,
             data=auc.wide)+
  coord_equal()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(experiment ~ test.samples)

auc.stats <- auc.dt[, list(
  median=median(auc),
  lo=quantile(auc, 0.25),
  hi=quantile(auc, 0.75)
  ), by=list(model, experiment, test.samples)]

