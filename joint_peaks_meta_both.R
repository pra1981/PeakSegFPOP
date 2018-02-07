library(ggplot2)
library(data.table)
sample.dir <- "labels/H3K4me3_TDH_immune/samples/monocyte/McGill0001"
jpeaks <- fread(file.path(sample.dir, "joint_peaks.bedGraph"))
coverage.bigWig <- file.path(sample.dir, "coverage.bigWig")
setnames(jpeaks, c("chrom", "peakStart", "peakEnd", "mean"))
library(PeakSegPipeline)

jpeaks[, peakBases := peakEnd-peakStart]
jpeaks[, log10peakBases := log10(peakBases)]
jpeaks[, log10.dist := log10peakBases-mean(log10peakBases)]
jpeaks[order(abs(log10.dist))]

n <- 5
list.dt <- jpeaks[, {
  zoomBases <- as.integer((peakEnd-peakStart)/2)
  from <- peakStart-zoomBases
  to <- peakEnd+zoomBases
  zoom.coverage <- readBigWig(coverage.bigWig, chrom, from, to)
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
write.table(mean.mat, file.path(sample.dir, "joint_peaks_meta.bedGraph"), col.names=FALSE, quote=FALSE)

mean.dt <- fread(file.path(sample.dir, "joint_peaks_meta.bedGraph"))
mean.mat <- as.matrix(mean.dt[,-1,with=FALSE])
rownames(mean.mat) <- mean.dt$V1

pc.fit <- princomp(mean.mat)
pc.dt <- data.table(
  pc.fit$scores[, 1:2],
  peakBases=list.dt[, peakEnd-peakStart])

ggplot()+
  geom_point(aes(Comp.1, Comp.2, color=log10(peakBases)), data=pc.dt)+
  scale_color_gradient(low="white", high="black")

ggplot()+
  geom_line(aes(box.i, mean, group=peakEnd), data=ten)+
  scale_x_continuous(
    "position relative to peak")

pattern <- paste0(
  "chr",
  "(?<chrom>[^:]")
i.list <- list(
  all=TRUE,
  iqr=list.dt[, which(899 < peakBases & peakBases < 2136)])

colnames(mean.mat) <- NULL
mean.tall <- data.table(melt(mean.mat))
setnames(mean.tall, c("peak", "pos", "mean"))
setkey(mean.tall, peak)

metapeak.dt.list <- list()
for(i.name in names(i.list)){
  i.vec <- i.list[[i.name]]
  name.vec <- rownames(mean.mat)[i.vec]
  i.stats <- mean.tall[name.vec, list(
    median=median(mean),
    q25=quantile(mean, 0.25),
    q75=quantile(mean, 0.75)
    ), by=list(pos)]
  metapeak.dt.list[[i.name]] <- data.table(i.name, i.stats)
}
metapeak.dt <- do.call(rbind, metapeak.dt.list)

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(i.name ~ .)+
  geom_ribbon(aes(pos, ymin=q25, ymax=q75), alpha=0.5, data=metapeak.dt)+
  geom_line(aes(pos, median), data=metapeak.dt)+
  scale_x_continuous(
    "position relative to peak",
    breaks=c(5.5, 15.5),
    labels=c("peakStart", "peakEnd"))
