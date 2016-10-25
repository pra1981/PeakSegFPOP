source("test_functions.R")

## Download trackDb.txt from github.
bigWig.part.vec <- c(
  "Input/MS010302",
  "bcell/MS010302",
  ## "Input/MS026601",
  ## "bcell/MS026601",
  ## "Input/MS002201",
  ## "kidney/MS002201",
  "Input/MS002202",
  "kidney/MS002202")
set.dir <- file.path("test", "input")
repos.url <- "https://raw.githubusercontent.com/tdhock/input-test-data/master/"
for(bigWig.part in bigWig.part.vec){
  bigWig.file <- file.path(set.dir, "samples", bigWig.part, "coverage.bigWig")
  bigWig.url <- paste0(repos.url, bigWig.part, ".bigwig")
  download.to(bigWig.url, bigWig.file)
}
labels.url <- paste0(repos.url, "kidney_bcell_Input_labels.txt")
labels.file <- file.path(set.dir, "labels", "kidney_bcell_Input_labels.txt")
download.to(labels.url, labels.file)
problems.bed <- file.path(set.dir, "problems.bed")
unlink(problems.bed)
system(paste("grep chr10 hg19_problems.bed >", problems.bed))

## Whole pipeline.
system(paste("bigWigToBedGraph", bigWig.file, "/dev/stdout|head"))
convert.cmd <- paste("Rscript pipeline.R", set.dir)
status <- system(convert.cmd)
test_that("pipeline script succeeds", {
  expect_equal(status, 0)
})
test_that("index.html is created", {
  expect_true(file.exists(file.path(set.dir, "index.html")))
})

## Plot peak predictions on input data set.
if(FALSE){

  joint.problems.dt <- fread(paste("cat", file.path(
    set.dir, "problems", "*", "jointProblems.bed")))
  setnames(joint.problems.dt, c("chrom", "problemStart", "problemEnd"))
  peaks.glob <- file.path(
    set.dir, "problems", "*", "jointProblems", "*", "peaks.bed")
  ##Sys.glob(peaks.glob)
  joint.peaks.dt <- fread(paste("cat", peaks.glob))
  setnames(
    joint.peaks.dt,
    c("chrom", "peakStart", "peakEnd", "sample.path", "mean"))
  joint.peaks.dt[, sample.id := sub(".*/", "", sample.path)]
  joint.peaks.dt[, sample.group := sub("/.*", "", sample.path)]
  chunks.dt <- fread(paste("cat", file.path(
    set.dir, "problems", "*", "chunks", "*", "chunk.bed")))
  setnames(chunks.dt, c("chrom", "chunkStart", "chunkEnd"))
  all.problems <- fread(file.path(set.dir, "problems.bed"))
  setnames(all.problems, c("chrom", "problemStart", "problemEnd"))
  data.start <- min(chunks.dt$chunkStart)
  data.end <- max(chunks.dt$chunkEnd)
  two.problems <- all.problems[!(
    problemEnd < data.start |
      data.end < problemStart),]
  both.problems <- rbind(
    data.table(two.problems, y="separate problems"),
    data.table(joint.problems.dt, y="joint problems"))

  limits.list <- list(
    c(34, 35),
    c(35, 36),
    c(36, 37),
    c(37, 38),
    c(38, 39))
  limits.dt <- data.table(matrix(unlist(limits.list)*1e6, ncol=2, byrow=TRUE))
  setnames(limits.dt, c("plotStart", "plotEnd"))
  limits.dt$y <- "plots in un-labeled regions"

  labels.bed.vec <- Sys.glob(file.path(
    set.dir, "samples", "*", "*", "labels.bed"))
  all.labels.list <- list()
  for(labels.bed in labels.bed.vec){
    sample.dir <- dirname(labels.bed)
    sample.id <- basename(sample.dir)
    group.dir <- dirname(sample.dir)
    sample.group <- basename(group.dir)
    sample.labels <- fread(labels.bed)
    setnames(sample.labels, c("chrom", "labelStart", "labelEnd", "annotation"))
    all.labels.list[[labels.bed]] <- data.table(
      sample.id, sample.group, sample.labels)
  }
  all.labels <- do.call(rbind, all.labels.list)

  input.labels <- all.labels[sample.group=="Input", list(
    prop.noPeaks=mean(annotation=="noPeaks")
  ), by=.(labelStart, labelEnd)]
  setkey(input.labels, labelStart, labelEnd)
  input.pred <- joint.peaks.dt[, list(
    n.Input=sum(sample.group=="Input")
  ), by=.(peakStart, peakEnd)]
  setkey(input.pred, peakStart, peakEnd)
  labeled.input <- foverlaps(input.pred, input.labels, nomatch=0L)
  thresh.dt <- labeled.input[, data.table(WeightedROC(
    n.Input, ifelse(prop.noPeaks==0, 1, -1)))]
  thresh.best <- thresh.dt[which.min(FP+FN),]
  ## threshold is smallest n.Input that is classified as non-specific.
  setkey(joint.peaks.dt, peakStart, peakEnd)
  peaks.with.counts <- input.pred[joint.peaks.dt]
  peaks.with.counts[, specificity := ifelse(
    n.Input >= thresh.best$threshold, "non-specific", "specific")]

  gg <- ggplot()+
    coord_cartesian(xlim=c(data.start, data.end)/1e3, expand=TRUE)+
    geom_segment(aes(
      problemStart/1e3, y,
      xend=problemEnd/1e3, yend=y),
      color="blue",
      data=both.problems)+
    geom_point(aes(
      problemStart/1e3, y),
      color="blue",
      data=both.problems)+
    geom_tallrect(aes(
      xmin=chunkStart/1e3, xmax=chunkEnd/1e3),
      alpha=0.1,
      data=chunks.dt)+
    geom_segment(aes(
      plotStart/1e3, y,
      xend=plotEnd/1e3, yend=y),
      data=limits.dt)+
    geom_point(aes(
      plotStart/1e3, y),
      data=limits.dt)+
    geom_point(aes(
      peakStart/1e3, sample.path, color=specificity),
      data=peaks.with.counts)+
    geom_segment(aes(
      peakStart/1e3, sample.path,
      xend=peakEnd/1e3, yend=sample.path, color=specificity),
      data=peaks.with.counts)+
    ylab("")+
    xlab("position on chr10 (kb = kilo bases)")+
    scale_color_manual(values=c("non-specific"="red", specific="deepskyblue"))
  png("figure-input-overview.png", res=100, h=200, w=1000)
  print(gg)
  dev.off()

  coverage.bigWig.vec <- Sys.glob(file.path(
    set.dir, "samples", "*", "*", "coverage.bigWig"))

  ann.colors <-
    c(noPeaks="#f6f4bf",
      peakStart="#ffafaf",
      peakEnd="#ff4c4c",
      peaks="#a445ee")

  for(limit.i in seq_along(limits.list)){
    limit.vec <- limits.list[[limit.i]]*1e6
    print(limit.vec)
    coverage.list <- list()
    for(coverage.bigWig in coverage.bigWig.vec){
      sample.dir <- dirname(coverage.bigWig)
      sample.id <- basename(sample.dir)
      group.dir <- dirname(sample.dir)
      sample.group <- basename(group.dir)
      cmd <- sprintf(
        "bigWigToBedGraph -chrom=chr10 -start=%d -end=%d %s /dev/stdout",
        limit.vec[1], limit.vec[2],
        coverage.bigWig)
      sample.coverage <- fread(cmd)
      setnames(sample.coverage, c("chrom", "chromStart", "chromEnd", "count"))
      coverage.list[[coverage.bigWig]] <- data.table(
        sample.id, sample.group, sample.coverage)
    }
    coverage <- do.call(rbind, coverage.list)

    show.peaks <- joint.peaks.dt[
      !(peakEnd < limit.vec[1] | limit.vec[2] < peakStart),]
    show.labels <- all.labels[
      !(labelEnd < limit.vec[1] | limit.vec[2] < labelStart),]

    gg <- ggplot()+
      theme_bw()+
      theme(panel.margin=grid::unit(0, "lines"))+
      facet_grid(sample.group + sample.id ~ ., scales="free")+
      geom_tallrect(aes(
        xmin=labelStart/1e3, 
        xmax=labelEnd/1e3,
        fill=annotation), 
        alpha=0.5,
        data=show.labels)+
      scale_fill_manual(values=ann.colors)+
      scale_y_continuous(
        "aligned read coverage",
        breaks=function(limits){
          lim <- floor(limits[2])
          if(lim==0){
            Inf
          }else{
            lim
          }
        })+
      scale_x_continuous(paste(
        "position on",
        coverage$chrom[1],
        "(kb = kilo bases)"))+
      geom_step(aes(
        chromStart/1e3, count),
        color="grey50",
        data=coverage)+
      coord_cartesian(xlim=limit.vec/1e3, expand=FALSE)
    if(nrow(show.peaks)){
      gg <- gg+
        geom_point(aes(
          peakStart/1e3, 0),
          color="deepskyblue",
          size=3,
          data=show.peaks)+
        geom_segment(aes(
          peakStart/1e3, 0,
          xend=peakEnd/1e3, yend=0),
          color="deepskyblue",
          size=3,
          data=show.peaks)
    }
    f <- sprintf("figure-input-%d-%d.png", limit.vec[1]/1e6, limit.vec[2]/1e6)
    png(f, res=100, width=1000, height=600)
    print(gg)
    dev.off()
    
  }

}  

