library(data.table)
library(flexclust)
library(ggplot2)
macs.default.peaks <- fread("macs.default.peaks.csv")
macs.default.peaks[, peakBases := peakEnd - peakStart]
macs.chrom.counts <- dcast(macs.default.peaks, experiment + sample.id ~ chrom, length)
setkey(macs.chrom.counts, experiment)

## Plot peak height and size for one sample...
peaks1 <- data.table(joint_peaks.bedGraph=Sys.glob("labels/H*_TDH_*/samples/monocyte/McGill0001/joint_peaks.bedGraph"))[, {
  peaks <- fread(joint_peaks.bedGraph)
  setnames(peaks, c("chrom", "peakStart", "peakEnd", "mean"))
  sample.dir <- dirname(joint_peaks.bedGraph)
  sample.id <- basename(sample.dir)
  group.dir <- dirname(sample.dir)
  samples.dir <- dirname(group.dir)
  proj.dir <- dirname(samples.dir)
  set.name <- basename(proj.dir)
  experiment <- sub("_.*", "", set.name)
  train.samples <- sub(".*_", "", set.name)
  data.table(experiment, train.samples, peaks)
}, by=list(joint_peaks.bedGraph)]
peaks1[, peakBases := peakEnd - peakStart]

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(experiment ~ train.samples)+
  geom_hex(aes(log10(peakBases), log10(mean)), data=peaks1)
    

## count peaks per chrom, compare clusterings.
for(experiment in c("H3K4me3", "H3K36me3")){
  macs.exp.dt <- macs.chrom.counts[experiment]
  setkey(macs.exp.dt, sample.id)
  peakseg.dir <- file.path("labels", paste0(experiment, "_TDH_other"))
  peaks.matrix.tsv <- file.path(peakseg.dir, "peaks_matrix.tsv")
  peakseg.dt <- fread(peaks.matrix.tsv)
  summary.tsv <- file.path(peakseg.dir, "peaks_summary.tsv")
  summary.dt <- fread(summary.tsv)
  peakseg.dt[, chrom := summary.dt$chrom]
  peakseg.tall.dt <- peakseg.dt[, {
    count.mat <- as.matrix(.SD[, !names(.SD) %in% c("peak", "chrom"), with=FALSE])
    peaks <- colSums(count.mat)
    sample.id <- sub(".*/", "", names(peaks))
    data.table(sample.path=names(peaks), peaks)
  }, by=list(chrom=factorChrom(chrom))]
  peakseg.wide.dt <- dcast(peakseg.tall.dt, sample.path ~ chrom, value.var="peaks")
  peakseg.wide.mat <- as.matrix(peakseg.wide.dt[, -1, with=FALSE])
  rownames(peakseg.wide.mat) <- peakseg.wide.dt$sample.path
  id2path <- peakseg.wide.dt$sample.path
  names(id2path) <- sub(".*/", "", id2path)
  macs.exp.mat <- as.matrix(macs.exp.dt[names(id2path), colnames(peakseg.wide.mat), with=FALSE])
  rownames(macs.exp.mat) <- id2path
  mat.list <- list(
    macs=macs.exp.mat,
    peakseg=peakseg.wide.mat)
  for(model.name in names(mat.list)){
    mat <- mat.list[[model.name]]
    d <- dist(mat, "manhattan")
    tree <- hclust(d, "average")
    plot(tree)
    guess.mat <- cutree(tree, h=tree$height)
    lab.vec <- factor(sub("/.*", "", rownames(guess.mat)))
    rand.vec <- apply(guess.mat, 2, flexclust::randIndex, lab.vec)
    plot(rand.vec)
    max(rand.vec)
    pc.fit <- princomp(mat)
    pc.dt <- with(pc.fit, data.table(scores[, c("Comp.1", "Comp.2")], label=lab.vec))
    ggplot()+
      geom_point(aes(Comp.1, Comp.2, color=label), data=pc.dt)
  }
}

## check peak size.
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(experiment ~ .)+
  geom_histogram(aes(log10(peakBases), ..density..), data=macs.default.peaks)

peakseg.dt <- data.table(summary.tsv=Sys.glob("labels/H*_TDH_*/peaks_summary.tsv"))[, {
  summary.dt <- fread(summary.tsv)
  data.dir <- basename(dirname(summary.tsv))
  experiment <- sub("_.*", "", data.dir)
  model <- paste0("PeakSeg.", sub(".*_", "", data.dir))
  summary.dt[, data.table(experiment, model, peakBases)]
}, by=list(summary.tsv)]

bases.dt <- rbind(
  data.table(peakseg.dt[, c("model", "experiment", "peakBases"), with=FALSE]),
  data.table(macs.default.peaks[, c("model", "experiment", "peakBases"), with=FALSE]))

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(experiment ~ .)+
  geom_freqpoly(aes(log10(peakBases), ..density.., group=model, color=model), data=bases.dt)

bases.stats <- bases.dt[, list(
  q25=quantile(peakBases, 0.25),
  median=as.numeric(median(peakBases)),
  q75=quantile(peakBases, 0.75)
  ), by=list(model, experiment)]

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(experiment ~ .)+
  geom_point(aes(median, model), data=bases.stats)+
  geom_segment(aes(
    q25, model,
    xend=q75, yend=model), data=bases.stats)+
  scale_x_log10("peak size (bases)")

setkey(bases.stats, experiment, model)
gg <- ggplot()+
  geom_point(aes(
    median, model, color=experiment),
             size=4,
             shape=1,
             data=bases.stats)+
  geom_segment(aes(
    q25, model,
    color=experiment,
    size=experiment,
    xend=q75, yend=model), data=bases.stats)+
  scale_x_log10("peak size (bases), median and quartiles")+
  scale_size_manual(values=c(H3K4me3=1, H3K36me3=2))
png("labels/H3K4me3_TDH_immune/figure-PeakSeg-macs-peak-size-comparison.png", 500, 150, res=100)
print(gg)
dev.off()

## test error.
setkey(macs.default.peaks, experiment, sample.id, model, chrom, peakStart, peakEnd)
error.dt <- data.table(labels.bed=Sys.glob("labels/H*_TDH_*/samples/*/*/labels.bed"))[, {
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
  joint_peaks.bedGraph <- file.path(
    "labels", train.proj.dir, "samples",
    basename(group.dir), sample.id,
    "joint_peaks.bedGraph")
  if(file.exists(joint_peaks.bedGraph)){
    PeakSeg.dt <- fread(joint_peaks.bedGraph)
    setnames(PeakSeg.dt, c("chrom", "chromStart", "chromEnd", "mean"))
    select.dt <- data.table(experiment, sample.id)
    pred.dt <- rbind(
      PeakSeg.dt[, data.table(model="PeakSeg", chrom, chromStart, chromEnd)],
      macs.default.peaks[select.dt, list(
        model, chrom, chromStart=peakStart, chromEnd=peakEnd)])
    browser()
    error.dt <- pred.dt[, {
      PeakError::PeakError(.SD, labels.dt)
    }, by=list(model)]
    data.table(experiment, test.samples, error.dt)
  }
}, by=list(labels.bed)]

error.totals <- error.dt[, {
  fp <- sum(fp)
  fn <- sum(fn)
  errors <- fp+fn
  data.table(
    fp, fn, errors,
    labels=.N,
    percent.errors=errors/.N,
    possible.fp=sum(possible.fp),
    possible.fn=sum(possible.tp))
}, by=list(experiment, test.samples, model)]
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(experiment ~ .)+
  geom_point(aes(
    percent.errors, model, color=test.samples),
             data=error.totals)

