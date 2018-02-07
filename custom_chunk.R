chunk.dir <- "labels/H3K36me3_TDH_immune/problems/chr10:18024675-38818835/chunks/chr10:18761902-22380580"
chunks.dir <- dirname(chunk.dir)
prob.dir <- dirname(chunks.dir)
prob.name <- basename(prob.dir)
probs.dir <- dirname(prob.dir)
proj.dir <- dirname(probs.dir)
problem.dir.vec <- Sys.glob(file.path(
  proj.dir, "samples", "*", "*", "problems", prob.name))
chunk <- fread(file.path(chunk.dir, "chunk.bed"))
setnames(chunk, c("chrom", "chunkStart", "chunkEnd"))
chunk[, chunkStart1 := chunkStart + 1L]
setkey(chunk, chunkStart1, chunkEnd)
labels <- fread(file.path(chunk.dir, "labels.tsv"))
cat("Read",
    nrow(labels),
    "labels.\n")
jointProblems <- fread(file.path(prob.dir, "jointProblems.bed"))
setnames(jointProblems, c("chrom", "problemStart", "problemEnd"))
jointProblems[, problemStart1 := problemStart + 1L]
jointProblems[, problem.name := sprintf(
                               "%s:%d-%d", chrom, problemStart, problemEnd)]
setkey(jointProblems, problemStart1, problemEnd)
probs.in.chunk <- foverlaps(jointProblems, chunk, nomatch=0L)
probs.in.chunk$sample.group <- "problems"
probs.in.chunk$sample.id <- "joint"
cat("Read",
    nrow(jointProblems),
    "joint problems, plotting",
    nrow(probs.in.chunk),
    "in chunk.\n")
grep.str <- "0007|0013|0022|0002"
problem.dir.vec <- grep(grep.str, problem.dir.vec, value=TRUE)
coverage.list <- list()
separate.peaks.list <- list()
for(sample.i in seq_along(problem.dir.vec)){
  problem.dir <- problem.dir.vec[[sample.i]]
  problems.dir <- dirname(problem.dir)
  sample.dir <- dirname(problems.dir)
  sample.id <- basename(sample.dir)
  group.dir <- dirname(sample.dir)
  sample.group <- basename(group.dir)
  sample.coverage <- fread(file.path(problem.dir, "coverage.bedGraph"))
  setnames(sample.coverage, c("chrom", "chromStart", "chromEnd", "count"))
  sample.peaks <- fread(file.path(problem.dir, "peaks.bed"))
  setnames(sample.peaks, c("chrom", "peakStart", "peakEnd", "status", "mean"))
  sample.peaks[, peakStart1 := peakStart + 1L]
  sample.coverage[, chromStart1 := chromStart + 1L]
  setkey(sample.coverage, chromStart1, chromEnd)
  setkey(sample.peaks, peakStart1, peakEnd)
  chunk.cov <- foverlaps(sample.coverage, chunk, nomatch=0L)
  chunk.peaks <- foverlaps(sample.peaks, chunk, nomatch=0L)
  coverage.list[[problem.dir]] <- data.table(
    sample.id, sample.group, chunk.cov)
  if(nrow(chunk.peaks)){
    separate.peaks.list[[problem.dir]] <- data.table(
      sample.id, sample.group, chunk.peaks)
  }
}
coverage <- do.call(rbind, coverage.list)
separate.peaks <- do.call(rbind, separate.peaks.list)
separate.peaks$peak.type <- "separate"
cat("Read",
    length(coverage.list),
    "samples of coverage.\n")
cat("Read",
    length(separate.peaks.list),
    "samples of separate peak predictions.\n")
joint.peaks.list <- list()
for(joint.i in 1:nrow(probs.in.chunk)){
  prob <- probs.in.chunk[joint.i,]
  tryCatch({
    peaks <- fread(file.path(
      prob.dir, "jointProblems", prob$problem.name, "peaks.bed"))
    setnames(peaks, c("chrom", "peakStart", "peakEnd", "sample.path", "mean"))
    peaks[, sample.id := sub(".*/", "", sample.path)]
    peaks[, sample.group := sub("/.*", "", sample.path)]
    joint.peaks.list[[prob$problem.name]] <- peaks[grepl(grep.str, sample.id),]
  }, error=function(e){
    NULL
  })
}
cat("Read",
    length(joint.peaks.list),
    "joint peak predictions.\n")
joint.peaks <- do.call(rbind, joint.peaks.list)
ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(sample.id + sample.group ~ ., scales="free")+
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
  ## geom_tallrect(aes(
  ##   xmin=problemStart/1e3,
  ##   xmax=problemEnd/1e3),
  ##   alpha=0.5,
  ##   size=3,
  ##   color="black",
  ##   fill=NA,
  ##   data=probs.in.chunk)+
  geom_tallrect(aes(
    xmin=chromStart/1e3, 
    xmax=chromEnd/1e3,
    fill=annotation), 
                alpha=0.5,
                data=labels[grepl(grep.str, sample.id),])+
  scale_fill_manual("label", values=ann.colors)+
  scale_color_manual(values=c(separate="black", joint="deepskyblue"))+
  scale_size_manual(values=c(separate=2, joint=3))+
  geom_step(aes(
    chromStart/1e3, count),
            data=coverage,
            color="grey50")
if(length(joint.peaks)){
  joint.peaks$peak.type <- "joint"
  gg <- gg+
    geom_point(aes(
      peakStart/1e3, 0,
      color=peak.type,
      size=peak.type),
               data=joint.peaks)+
                 geom_segment(aes(
                   peakStart/1e3, 0,
                   xend=peakEnd/1e3, yend=0,
                   color=peak.type,
                   size=peak.type),
                              data=joint.peaks)
}
gg <- gg+
  geom_segment(aes(
    peakStart/1e3, 0,
    xend=peakEnd/1e3, yend=0,
    color=peak.type,
    size=peak.type),
               data=separate.peaks)+
  geom_point(aes(
    peakStart/1e3, 0,
    color=peak.type,
    size=peak.type),
             data=separate.peaks)

n.rows <- length(coverage.list) + 2
mypng <- function(base, g){
  f <- file.path(chunk.dir, base)
  cat("Writing ",
      f,
      "\n", sep="")
  png(f, res=100, width=1000, height=100*n.rows)
  print(g)
  dev.off()
  thumb.png <- sub(".png$", "-thumb.png", f)
  cmd <- sprintf("convert %s -resize 230 %s", f, thumb.png)
  system(cmd)
}
gg.zoom <- gg+
  coord_cartesian(
    xlim=chunk[, c(chunkStart, chunkEnd)/1e3],
    expand=FALSE)
mypng("figure-predictions-custom.png", gg.zoom)
