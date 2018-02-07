library(ggplot2)
library(data.table)
library(PeakSegPipeline)

sample.dir <- "labels/H3K4me3_TDH_immune/samples/monocyte/McGill0001"
sample.dir <- "labels/H3K36me3_TDH_immune/samples/monocyte/McGill0001"

source("getMeanMat.R")
## This is for creating macs2.broad_peaks.bedGraph files, for input to
## getMeanMat above.
macs.default.peaks <- fread("macs.default.peaks.csv")
macs.default.peaks[sample.id=="McGill0001", {
  set.name <- paste0(experiment, "_TDH_immune")
  sample.dir <- file.path("labels", set.name, "samples/monocyte/McGill0001")
  peaks.bedGraph <- file.path(sample.dir, paste0(model, "_peaks.bedGraph"))
  fwrite(.SD[, list(
    chrom, peakStart, peakEnd, mean=NA)],
         peaks.bedGraph,
         sep="\t",
         col.names=FALSE)
}, by=list(model, experiment)]

pc.fit <- princomp(mean.mat)
pc.dt <- data.table(
  pc.fit$scores[, 1:2],
  peakBases=k36.dt[, peakEnd-peakStart])

ggplot()+
  geom_point(aes(Comp.1, Comp.2, color=log10(peakBases)), data=pc.dt)+
  scale_color_gradient(low="white", high="black")

ggplot()+
  geom_line(aes(box.i, mean, group=peakEnd), data=ten)+
  scale_x_continuous(
    "position relative to peak")

pattern <- paste0(
  "(?<chrom>chr[^:]+)",
  ":",
  "(?<peakStart>[^-]+)",
  "-",
  "(?<peakEnd>.*)")
size.dt <- data.table(str_match_named(rownames(mean.mat), pattern, list(
  peakStart=as.integer,
  peakEnd=as.integer)))
size.dt[, peakBases := peakEnd - peakStart]
quantile.vec <- quantile(size.dt$peakBases)
i.list <- list(
  all=TRUE,
  iqr=size.dt[, which(quantile.vec[["25%"]] < peakBases & peakBases < quantile.vec[["75%"]])])

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


exp.stats.dt <- data.table(joint_peaks_mat.csv=Sys.glob("labels/H*_TDH_immune/samples/monocyte/McGill0001/joint_peaks_mat.csv"))[, {
  mean.mat <- readMeanMat(joint_peaks_mat.csv)
  mean.tall <- data.table(melt(mean.mat))
  setnames(mean.tall, c("peak", "pos", "mean"))
  i.stats <- mean.tall[, list(
    median=median(mean),
    q25=quantile(mean, 0.25),
    q75=quantile(mean, 0.75)
    ), by=list(pos)]
  sample.dir <- dirname(joint_peaks_mat.csv)
  group.dir <- dirname(sample.dir)
  samples.dir <- dirname(group.dir)
  set.dir <- dirname(samples.dir)
  set.name <- basename(set.dir)
  experiment <- sub("_.*", "", set.name)
  data.table(experiment, i.stats)
}, by=list(joint_peaks_mat.csv)]

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(experiment ~ ., scales="free")+
  geom_ribbon(aes(pos, ymin=q25, ymax=q75), alpha=0.5, data=exp.stats.dt)+
  geom_line(aes(pos, median), data=exp.stats.dt)+
  scale_x_continuous(
    "position relative to peak",
    breaks=c(5.5, 15.5),
    labels=c("peakStart", "peakEnd"))+
  scale_y_continuous("aligned read coverage (median and quartiles)")

## Parallel version.
peaks.bedGraph.vec <- Sys.glob("labels/H*_TDH_immune/samples/monocyte/McGill0001/*_peaks.bedGraph")
peaks.bedGraph.vec <- file.path(
  "labels/H3K27ac-H3K4me3_TDHAM_BP/samples",
  c("GR4_H3K27ac/S00G5U_WTSI",
    "GR12_H3K27ac/S00XGE_WTSI",
    "GR12_H3K27ac/S00XHC_WTSI"),
  "joint_peaks.bedGraph")
library(parallel)
options(mc.cores=4)
exp.stats.dt.list <- mclapply(peaks.bedGraph.vec, function(peaks.bedGraph){
  print(peaks.bedGraph)
  mean.mat <- getMeanMat(peaks.bedGraph)
  mean.tall <- data.table(melt(mean.mat))
  setnames(mean.tall, c("peak", "pos", "mean"))
  i.stats <- mean.tall[, list(
    median=median(mean),
    q25=quantile(mean, 0.25),
    q75=quantile(mean, 0.75),
    count=.N
    ), by=list(pos)]
  sample.dir <- dirname(peaks.bedGraph)
  group.dir <- dirname(sample.dir)
  samples.dir <- dirname(group.dir)
  set.dir <- dirname(samples.dir)
  set.name <- basename(set.dir)
  experiment <- sub("_.*", "", set.name)
  model <- sub("_.*", "", basename(peaks.bedGraph))
  data.table(experiment, model, i.stats)
})
exp.stats.dt <- do.call(rbind, exp.stats.dt.list)

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(experiment ~ model, scales="free")+
  geom_ribbon(aes(pos, ymin=q25, ymax=q75), alpha=0.5, data=exp.stats.dt)+
  geom_line(aes(pos, median), data=exp.stats.dt)+
  scale_x_continuous(
    "position relative to peak",
    breaks=c(5.5, 15.5),
    labels=c("peakStart", "peakEnd"))+
  scale_y_continuous("aligned read coverage (median and quartiles)")+
  geom_text(aes(
    5.5, -Inf, label=paste(count, "peaks")), data=exp.stats.dt[pos==5],
            hjust=0,
            vjust=-0.5)
print(gg)

png("figure-joint-peaks-meta.png", 800, 500, res=100)
print(gg)
dev.off()

