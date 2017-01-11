works_with_R("3.3.2",
             "tdhock/namedCapture@05175927a45c301a18e8c6ebae67ea39a842d264")
trackDb.lines.vec <- readLines("trackDb_2epigenomes.txt")

pattern <- paste0(
  "bigDataUrl",
  " ",
  "(?<suffix>.*[.]bigWig)",
  "\n")
trackDb.txt <- paste(trackDb.lines.vec, collapse="\n")
match.mat <- str_match_all_named(trackDb.txt, pattern)[[1]]
prefix <- "http://epigenomesportal.ca"
suffix.vec <- match.mat[, "suffix"]
dir.create("trackDb_2epigenomes", showWarnings=FALSE)
for(suffix in suffix.vec){
  suffix.base <- basename(suffix)
  local.path <- file.path("trackDb_2epigenomes", suffix.base)
  if(!file.exists(local.path)){
    u <- paste0(prefix, suffix)
    download.file(u, local.path)
  }
}

hg19.probs <- fread("hg19_problems.bed")
setnames(hg19.probs, c("chrom", "chromStart", "chromEnd"))
hg19.probs[, chrom.fac := factorChrom(chrom)]
chrom.levs <- levels(hg19.probs$chrom.fac)

bigWig.vec <- Sys.glob(file.path("trackDb_2epigenomes", "*.bigWig"))
for(bigWig in bigWig.vec){
  bedGraph <- sub("bigWig$", "bedGraph", bigWig)
  top.file <- sub("unstranded", "q0.999", bedGraph)
  if(!file.exists(top.file)){
    if(!file.exists(bedGraph)){
      cmd <- paste("bigWigToBedGraph", bigWig, bedGraph)
      print(cmd)
      system(cmd)
    }
    print(bedGraph)
    coverage.dt <- fread(bedGraph)
    gc()
    setnames(coverage.dt, c("chrom", "chromStart", "chromEnd", "coverage.norm"))
    min.nonzero <- coverage.dt[0 < coverage.norm, min(coverage.norm)]
    coverage.dt[, count := as.integer(coverage.norm / min.nonzero)]
    coverage.dt[, quantile(count, seq(0.99, 1, by=0.001))]
    large.count <- quantile(coverage.dt$count, 0.999)
    large.dt <- coverage.dt[large.count < count, ]
    large.dt[, chrom.fac := factor(chrom, chrom.levs)]
    large.dt[, bases.to.next := c(chromStart[-1] - chromEnd[-.N], NA)]
    gc()

    ggplot()+
      scale_y_discrete(limits=chrom.levs)+
      geom_segment(aes(chromStart/1e6, chrom.fac, xend=chromEnd/1e6, yend=chrom.fac),
                   color="grey",
                   size=1,
                   data=hg19.probs)+
      ## geom_point(aes(chromStart/1e6, chrom.fac, yend=chrom.fac),
      ##              color="grey50",
      ##              data=hg19.probs)+
      geom_point(aes(chromStart/1e6, chrom.fac), data=large.dt, shape=1)

    fwrite(large.dt[, .(chrom, chromStart, chromEnd, count)], top.file, sep="\t")
    
  }
}

large.file.vec <- Sys.glob("trackDb_2epigenomes/*q0.999.bedGraph")
pattern <- paste0(
  "trackDb_2epigenomes/",
  "[0-9]+",
  ".CEEHRC.",
  "(?<donor>[^.]+)",
  ".",
  "(?<experiment>[^.]+)",
  ".signal_q0.999.bedGraph")

large.list <- list()
for(large.file in large.file.vec){
  large.dt <- fread(large.file)
  meta <- str_match_named(large.file, pattern)
  if(any(is.na(meta))){
    stop("NA!")
  }
  large.list[[large.file]] <- data.table(meta, large.dt)
}
large <- do.call(rbind, large.list)

large[, type := ifelse(experiment=="Input", "Input", "experiment")]
library(PeakSegDP)
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(chrom ~ .)+
  geom_tallrect(aes(xmin=chromStart/1e6, xmax=chromEnd/1e6),
               color="grey",
               size=1,
               data=hg19.probs)+
  geom_point(aes(chromStart/1e6, type, color=donor), data=large, shape=1)
