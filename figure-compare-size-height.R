library(parallel)
library(ggplot2)
library(data.table)
library(PeakSegPipeline)
source("getMeanMat.R")
options(mc.cores=4)

## experiment 1: compare three models and two experiments for the same sample.
peaks.bedGraph.vec <- Sys.glob("labels/H*_TDH_immune/samples/monocyte/McGill0001/*_peaks.bedGraph")
exp.stats.dt.list <- mclapply(peaks.bedGraph.vec, function(peaks.bedGraph){
  print(peaks.bedGraph)
  mean.mat <- getMeanMat(peaks.bedGraph)
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
  size.dt[, peakHeight := rowMeans(mean.mat[, -c(1:5, 16:20)]) ]
  sample.dir <- dirname(peaks.bedGraph)
  group.dir <- dirname(sample.dir)
  samples.dir <- dirname(group.dir)
  set.dir <- dirname(samples.dir)
  set.name <- basename(set.dir)
  experiment <- sub("_.*", "", set.name)
  model <- sub("_.*", "", basename(peaks.bedGraph))
  data.table(experiment, model, size.dt)
})
exp.stats.dt <- do.call(rbind, exp.stats.dt.list)

mean.dt <- exp.stats.dt[, list(
  mean.peakBases=mean(peakBases),
  mean.peakHeight=mean(peakHeight),
  median.peakBases=median(peakBases),
  median.peakHeight=median(peakHeight)
  ), by=list(experiment, model)]

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(experiment ~ model)+
  geom_bin2d(aes(log10(peakBases), log10(peakHeight), fill=..density..), data=exp.stats.dt)+
  scale_fill_gradient(low="white", high="black")+
  ## geom_point(aes(
  ##   log10(mean.peakBases), log10(mean.peakHeight)
  ##   ), data=mean.dt,
  ##            color="red")+
  geom_point(aes(
    log10(median.peakBases), log10(median.peakHeight)
    ), data=mean.dt,
             shape=1,
             color="red")+
  geom_text(aes(
    4, 0,
    label=paste0(
      "medians:
height=", round(median.peakHeight, 1), " reads
width=", round(median.peakBases), " bases")),
            data=mean.dt, color="red")
##print(gg)

png("figure-compare-size-height.png", 800, 500, res=100)
print(gg)
dev.off()

## experiment 2: compare three samples of the same experiment type,
## with varying coverage.
peaks2.bedGraph.vec <- file.path(
  "labels/H3K27ac-H3K4me3_TDHAM_BP/samples",
  c("GR4_H3K27ac/S00G5U_WTSI",#0.05
    "GR12_H3K27ac/S00XGE_WTSI",#0.27
    "GR12_H3K27ac/S00YBK_WTSI"#1.17
    ##"GR12_H3K27ac/S00XHC_WTSI"#1.24 BAD?
    ),
  "joint_peaks.bedGraph")
cov.dt <- fread("labels/H3K27ac-H3K4me3_TDHAM_BP/coverage.csv", integer64="double")
setkey(cov.dt, input.bigWig)
exp.stats.dt.list <- lapply(peaks2.bedGraph.vec, function(peaks.bedGraph){
  print(peaks.bedGraph)
  mean.mat <- getMeanMat(peaks.bedGraph)
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
  size.dt[, peakHeight := rowMeans(mean.mat[, -c(1:5, 16:20)]) ]
  sample.dir <- dirname(peaks.bedGraph)
  group.dir <- dirname(sample.dir)
  samples.dir <- dirname(group.dir)
  set.dir <- dirname(samples.dir)
  set.name <- basename(set.dir)
  key <- file.path(sub(".*BP/", "", sample.dir), "coverage.bigWig")
  some.cov <- cov.dt[key, list(total.coverage, total.bases, mean.coverage)]
  mean.coverage <- round(some.cov$mean.coverage, 2)
  sample.id <- basename(sample.dir)
  experiment <- sub("_.*", "", set.name)
  model <- sub("_.*", "", basename(peaks.bedGraph))
  data.table(
    experiment, model, sample.id,
    ##some.cov,
    mean.coverage,
    size.dt)
})
exp.stats.dt <- do.call(rbind, exp.stats.dt.list)
exp.stats.dt[, mean.coverage.str := paste0(mean.coverage, "x")]
mean.dt <- exp.stats.dt[, list(
  mean.peakBases=mean(peakBases),
  mean.peakHeight=mean(peakHeight),
  median.peakBases=median(peakBases),
  median.peakHeight=median(peakHeight)
  ), by=list(mean.coverage, mean.coverage.str)]

gg <- ggplot()+
  ggtitle("cellType=GR, experiment=H3K27ac, center=WTSI")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ mean.coverage.str)+
  geom_bin2d(aes(log10(peakBases), log10(peakHeight), fill=..density..), data=exp.stats.dt)+
  scale_fill_gradient(low="white", high="black")+
  ## geom_point(aes(
  ##   log10(mean.peakBases), log10(mean.peakHeight)
  ##   ), data=mean.dt,
  ##            color="red")+
  geom_point(aes(
    log10(median.peakBases), log10(median.peakHeight)
    ), data=mean.dt,
             shape=1,
             color="red")+
  geom_text(aes(
    3.25, 1.5,
    label=paste0(
      "medians:
height=", round(median.peakHeight, 1), " reads
width=", round(median.peakBases), " bases")),
            vjust=1,
            data=mean.dt, color="red")+
  coord_cartesian(xlim=c(2.5, 4), ylim=c(-0.5, 1.5))
##print(gg)
png("figure-compare-size-height-GR-H3K27ac-WTSI.png", 800, 500, res=100)
print(gg)
dev.off()

## experiment 3: same thing but for other samples picked from
## scatterplot. (from NCMLS instead of WTSI)
peaks3.bedGraph.vec <- file.path(
  "labels/H3K27ac-H3K4me3_TDHAM_BP/samples",
  c("GR1_H3K27ac/S0026A_NCMLS",#0.05
    "GR2_H3K27ac/S003JH_NCMLS",#0.27
    "GR3_H3K27ac/S006XE_NCMLS"#1.17
    ##"GR12_H3K27ac/S00XHC_WTSI"#1.24 BAD?
    ),
  "joint_peaks.bedGraph")
cov.dt <- fread("labels/H3K27ac-H3K4me3_TDHAM_BP/coverage.csv", integer64="double")
setkey(cov.dt, input.bigWig)
exp.stats.dt.list <- lapply(peaks3.bedGraph.vec, function(peaks.bedGraph){
  print(peaks.bedGraph)
  mean.mat <- getMeanMat(peaks.bedGraph)
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
  size.dt[, peakHeight := rowMeans(mean.mat[, -c(1:5, 16:20)]) ]
  sample.dir <- dirname(peaks.bedGraph)
  group.dir <- dirname(sample.dir)
  samples.dir <- dirname(group.dir)
  set.dir <- dirname(samples.dir)
  set.name <- basename(set.dir)
  key <- file.path(sub(".*BP/", "", sample.dir), "coverage.bigWig")
  some.cov <- cov.dt[key, list(total.coverage, total.bases, mean.coverage)]
  mean.coverage <- round(some.cov$mean.coverage, 2)
  sample.id <- basename(sample.dir)
  experiment <- sub("_.*", "", set.name)
  model <- sub("_.*", "", basename(peaks.bedGraph))
  data.table(
    experiment, model, sample.id,
    ##some.cov,
    mean.coverage,
    size.dt)
})
exp.stats.dt <- do.call(rbind, exp.stats.dt.list)

exp.stats.dt[, mean.coverage.str := paste0(mean.coverage, "x")]
mean.dt <- exp.stats.dt[, list(
  mean.peakBases=mean(peakBases),
  mean.peakHeight=mean(peakHeight),
  median.peakBases=median(peakBases),
  median.peakHeight=median(peakHeight)
  ), by=list(mean.coverage, mean.coverage.str)]

gg <- ggplot()+
  ggtitle("cellType=GR, experiment=H3K27ac, center=NCMLS")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ mean.coverage.str)+
  geom_bin2d(aes(log10(peakBases), log10(peakHeight), fill=..density..), data=exp.stats.dt)+
  scale_fill_gradient(low="white", high="black")+
  ## geom_point(aes(
  ##   log10(mean.peakBases), log10(mean.peakHeight)
  ##   ), data=mean.dt,
  ##            color="red")+
  geom_point(aes(
    log10(median.peakBases), log10(median.peakHeight)
    ), data=mean.dt,
             shape=1,
             color="red")+
  geom_text(aes(
    3.25, -0.5,
    label=paste0(
      "medians:
height=", round(median.peakHeight, 1), " reads
width=", round(median.peakBases), " bases")),
            vjust=0,
            data=mean.dt, color="red")+
  coord_cartesian(ylim=c(-0.5, 1.5), xlim=c(2.5, 4))
##print(gg)
png("figure-compare-size-height-GR-H2K27ac-NCMLS.png", 800, 500, res=100)
print(gg)
dev.off()

## experiment 4: same thing but for other experiment. 
peaks4.bedGraph.vec <- file.path(
  "labels/H3K27ac-H3K4me3_TDHAM_BP/samples",
  c("GR3_H3K4me3/S00K4G_WTSI",#0.05
    "GR5_H3K4me3/S00PVG_WTSI",#0.27
    "GR1_H3K4me3/S00FPV_WTSI"#1.17
    ##"GR12_H3K27ac/S00XHC_WTSI"#1.24 BAD?
    ),
  "joint_peaks.bedGraph")
cov.dt <- fread("labels/H3K27ac-H3K4me3_TDHAM_BP/coverage.csv", integer64="double")
setkey(cov.dt, input.bigWig)
exp.stats.dt.list <- lapply(peaks4.bedGraph.vec, function(peaks.bedGraph){
  print(peaks.bedGraph)
  mean.mat <- getMeanMat(peaks.bedGraph)
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
  size.dt[, peakHeight := rowMeans(mean.mat[, -c(1:5, 16:20)]) ]
  sample.dir <- dirname(peaks.bedGraph)
  group.dir <- dirname(sample.dir)
  samples.dir <- dirname(group.dir)
  set.dir <- dirname(samples.dir)
  set.name <- basename(set.dir)
  key <- file.path(sub(".*BP/", "", sample.dir), "coverage.bigWig")
  some.cov <- cov.dt[key, list(total.coverage, total.bases, mean.coverage)]
  mean.coverage <- round(some.cov$mean.coverage, 2)
  sample.id <- basename(sample.dir)
  experiment <- sub("_.*", "", set.name)
  model <- sub("_.*", "", basename(peaks.bedGraph))
  data.table(
    experiment, model, sample.id,
    ##some.cov,
    mean.coverage,
    size.dt)
})
exp.stats.dt <- do.call(rbind, exp.stats.dt.list)
exp.stats.dt[, mean.coverage.str := paste0(mean.coverage, "x")]
mean.dt <- exp.stats.dt[, list(
  mean.peakBases=mean(peakBases),
  mean.peakHeight=mean(peakHeight),
  median.peakBases=median(peakBases),
  median.peakHeight=median(peakHeight)
  ), by=list(mean.coverage, mean.coverage.str)]

gg <- ggplot()+
  ggtitle("cellType=GR, experiment=H3K4me3, center=WTSI")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ mean.coverage.str)+
  geom_bin2d(aes(log10(peakBases), log10(peakHeight), fill=..density..), data=exp.stats.dt)+
  scale_fill_gradient(low="white", high="black")+
  ## geom_point(aes(
  ##   log10(mean.peakBases), log10(mean.peakHeight)
  ##   ), data=mean.dt,
  ##            color="red")+
  geom_point(aes(
    log10(median.peakBases), log10(median.peakHeight)
    ), data=mean.dt,
             shape=1,
             color="red")+
  coord_cartesian(ylim=c(-0.5, 1.5), xlim=c(2.5, 4))+
  geom_text(aes(
    3.25, -0.5,
    label=paste0(
      "medians:
height=", round(median.peakHeight, 1), " reads
width=", round(median.peakBases), " bases")),
            vjust=0,
            data=mean.dt, color="red")
##print(gg)

png("figure-compare-size-height-GR-H3K4me3-WTSI.png", 800, 500, res=100)
print(gg)
dev.off()


## experiment 5: scatterplot of all Blueprint H3K27ac and H3K4me3 samples.
peaks.bedGraph.vec <- normalizePath(Sys.glob(
  "~/PeakSegFPOP/labels/H3K27ac-H3K4me3_TDHAM_BP/samples/*/*/joint_peaks.bedGraph"))
already.done <- file.exists(sub("peaks.bedGraph", "peaks_mat.csv",  peaks.bedGraph.vec))
done.bedGraph.vec <- peaks.bedGraph.vec[already.done]
exp.stats.dt.list <- mclapply(seq_along(done.bedGraph.vec), function(peaks.i){
  peaks.bedGraph <- done.bedGraph.vec[[peaks.i]]
  cat(sprintf("%4d / %4d %s\n", peaks.i, length(done.bedGraph.vec), peaks.bedGraph))
  mean.mat <- getMeanMat(peaks.bedGraph)
  peak.pattern <- paste0(
    "(?<chrom>chr[^:]+)",
    ":",
    "(?<peakStart>[^-]+)",
    "-",
    "(?<peakEnd>.*)")
  size.dt <- data.table(str_match_named(rownames(mean.mat), peak.pattern, list(
    peakStart=as.integer,
    peakEnd=as.integer)))
  size.dt[, peakBases := peakEnd - peakStart]
  size.dt[, peakHeight := rowMeans(mean.mat[, -c(1:5, 16:20)]) ]
  sample.dir <- dirname(peaks.bedGraph)
  group.dir <- dirname(sample.dir)
  samples.dir <- dirname(group.dir)
  set.dir <- dirname(samples.dir)
  set.name <- basename(set.dir)
  sample.id <- basename(sample.dir)
  experiments <- sub("_.*", "", set.name)
  group.id <- basename(group.dir)
  experiment <- sub(".*_", "", group.id)
  model <- sub("_.*", "", basename(peaks.bedGraph))
  data.table(experiment, sample.id, group.id, size.dt)
})
exp.stats.dt <- do.call(rbind, exp.stats.dt.list)

sample.pattern <- paste0(
  "^",
  "(?<person>[^_]+)",
  "_",
  "(?<center>[^_]+)",
  "(?:_",
  "(?<cellType>.*)",
  ")?",
  "$")
sample.match.mat <- str_match_named(exp.stats.dt$sample.id, sample.pattern)
stopifnot(all(!is.na(sample.match.mat)))
exp.stats.dt[, center := sample.match.mat[, "center" ] ]

cov.dt[, sample.id := basename(dirname(input.bigWig))]
cov.dt[, group.id := basename(dirname(dirname(input.bigWig)))]
setkey(cov.dt, sample.id, group.id)

sample.height.dt <- exp.stats.dt[, list(
  median=median(peakHeight),
  q25=quantile(peakHeight, 0.25),
  q75=quantile(peakHeight, 0.75)
  ), by=list(sample.id,group.id,experiment,center)]
setkey(sample.height.dt, sample.id, group.id)

sample.stats.dt <- cov.dt[sample.height.dt]

plotted.dt <- data.table(
  joint_peaks.bedGraph=c(peaks2.bedGraph.vec, peaks3.bedGraph.vec, peaks4.bedGraph.vec),
  plot.i=rep(2:4, each=3))
plotted.dt[, sample.id := basename(dirname(joint_peaks.bedGraph))]
plotted.dt[, group.id := basename(dirname(dirname(joint_peaks.bedGraph)))]
setkey(plotted.dt, sample.id, group.id)
plotted.stats.dt <- sample.stats.dt[plotted.dt]

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ center)+
  geom_point(aes(mean.coverage, median, color=experiment), shape=1, data=sample.stats.dt)+
  scale_x_log10("mean sample coverage (aligned reads/base)")+
  scale_y_log10("median peak height (aligned reads/base)")
png("figure-compare-size-height-scatter.png", 800, 500, res=100)
print(gg)
dev.off()
gg <- gg+
  geom_point(aes(mean.coverage, median), size=1,data=plotted.stats.dt)+
  geom_line(aes(mean.coverage, median, group=plot.i), data=plotted.stats.dt)
png("figure-compare-size-height-scatter-emph.png", 800, 500, res=100)
print(gg)
dev.off()
  
