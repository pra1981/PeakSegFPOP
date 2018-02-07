library(data.table)
library(flexclust)
library(ggplot2)
set.vec <- c(
  "H3K4me3_TDH_immune",
  "H3K36me3_TDH_immune",
  "H3K4me1_TDH_BP",
  "H3K9me3_TDH_BP",
  "H3K27ac_TDH_some",
  "H3K27me3_TDH_some")

## experiment 1: plot distribution of peak sizes in two samples per data set.
bg.dt <- data.table(set.name=set.vec)[, {
  data.table(joint_peaks.bedGraph=Sys.glob(file.path("labels", set.name, "samples/*/*/joint_peaks.bedGraph")))
}, by=set.name]
some.bg <- bg.dt[, {
  data.table(experiment=sub("_.*", "", set.name), .SD[1:2])
}, by=set.name]
stats.dt <- some.bg[, {
  peaks <- fread(joint_peaks.bedGraph)
  setnames(peaks, c("chrom", "peakStart", "peakEnd", "mean"))
  peaks[, peakBases := peakEnd-peakStart]
  peaks[, list(
    q25=quantile(peakBases, 0.25),
    median=as.numeric(median(peakBases)),
    q75=quantile(peakBases, 0.75)
  )]
}, by=list(experiment, joint_peaks.bedGraph)]
stats.dt[, id := sub(".*samples/", "", sub("/joint_peaks.bedGraph", "", joint_peaks.bedGraph))]
exp.sizes <- stats.dt[, list(mean=mean(median)), by=experiment][order(mean)]
stats.dt[, exp.fac := factor(experiment, exp.sizes$experiment)]

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(exp.fac ~ ., scales="free", space="free")+
  geom_point(aes(median, id), data=stats.dt)+
  geom_segment(aes(
    q25, id,
    xend=q75, yend=id), data=stats.dt)+
  scale_x_log10("peak size (bases)")+
  ylab("sample ID")

## experiment 2: plot distribution of peak sizes in all joint peak calls.
summary.tsv.vec <- Sys.glob("labels/*/peaks_summary.tsv")
stats.dt <- data.table(summary.tsv=summary.tsv.vec)[, {
  summary.dt <- fread(summary.tsv)
  summary.dt[, list(
    set.name=basename(dirname(summary.tsv)),
    q25=quantile(peakBases, 0.25),
    median=as.numeric(median(peakBases)),
    q75=quantile(peakBases, 0.75)
  )]
}, by=summary.tsv]
set.sizes <- stats.dt[, list(mean=mean(median)), by=set.name][order(mean)]
stats.dt[, set.fac := factor(set.name, set.sizes$set.name)]
stats.dt[, samples := sub(".*_", "", set.name)]
stats.dt[, center := ifelse(
  samples=="BP", "BluePrint", ifelse(
  samples=="adipose", "unknown", ifelse(
  samples=="ENCODE", "ENCODE", "CEEHRC")))]

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_point(aes(median, set.fac, color=center), data=stats.dt)+
  geom_segment(aes(
    q25, set.fac,
    color=center,
    xend=q75, yend=set.fac), data=stats.dt)+
  scale_x_log10("peak size (bases)", limits=c(1e2, NA), breaks=10^seq(2, 5, by=1))+
  ylab("data set")
png("figure-compare-size-centers.png", 600, 200, res=100)
print(gg)
dev.off()
system("cp figure-compare-size-centers.png labels/H3K36me3_TDH_ENCODE")

  
## experiment 3: plot distribution of peak sizes of peaks present in
## all samples.
summary.tsv.vec <- Sys.glob("labels/*/peaks_summary.tsv")
some.vec <- grep("-", summary.tsv.vec, invert=TRUE, value=TRUE)
stats.dt <- data.table(summary.tsv=some.vec)[, {
  summary.dt <- fread(summary.tsv)
  no.input <- if("n.Input" %in% names(summary.dt)){
    summary.dt[n.Input==0]
  }else{
    summary.dt[, n.samples := n.samples.up]
    summary.dt[Input.up==FALSE]
  }
  max.samples <- max(no.input$n.samples)
  no.input[n.samples==max.samples, list(
    set.name=basename(dirname(summary.tsv)),
             max.samples,
    q25=quantile(peakBases, 0.25),
    median=as.numeric(median(peakBases)),
    q75=quantile(peakBases, 0.75),
             min=min(peakBases),
             max=max(peakBases),
               peaks=.N
  )]
}, by=summary.tsv]
set.sizes <- stats.dt[, list(mean=mean(median)), by=set.name][order(mean)]
stats.dt[, set.fac := factor(set.name, set.sizes$set.name)]
stats.dt[, samples := sub(".*_", "", set.name)]
stats.dt[, center := ifelse(samples=="BP", "BluePrint", ifelse(samples=="adipose", "unknown", "CEEHRC"))]
text.dt <- melt(stats.dt, measure.vars=c("q25", "q75"), id.vars=c("set.fac", "center"))

gg <- ggplot()+
  ggtitle("Inter-quartile range (25-75%) of sizes of peaks up in all samples")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_text(aes(
    value, set.fac, color=center,
    hjust=ifelse(variable=="q75", 0, 1),
    label=paste0(
      ifelse(variable=="q75", " ", ""),
      round(value),
      ifelse(variable=="q75", "", " "))),
    vjust=-0.2,
    data=text.dt)+
  geom_text(aes(
    median, set.fac, color=center,
    label=sprintf("N=%d peaks up in %d samples", as.integer(peaks/2), max.samples)),
            vjust=1.3,
            data=stats.dt)+
  geom_point(aes(
    median, set.fac, color=center),
             shape=1,
             data=stats.dt)+
  geom_segment(aes(
    q25, set.fac,
    color=center,
    xend=q75, yend=set.fac),
               size=1,
               data=stats.dt)+
  ## geom_segment(aes(
  ##   min, set.fac,
  ##   color=center,
  ##   xend=max, yend=set.fac),
  ##              size=0.5,
  ##              data=stats.dt)+
  scale_x_log10("peak size (bases)", limits=c(1e1, 1e6), breaks=10^seq(2, 5, by=1))+
  ylab("data set")
png("figure-compare-size-centers-allup.png", 1000, 500, res=100)
print(gg)
dev.off()
system("cp figure-compare-size-centers-allup.png labels/CTCF_TDH_ENCODE/")
  
