library(data.table)
problem.dir <- "labels/H3K4me3_TDH_immune/samples/tcell/McGill0025/problems/chr18:18521000-52059136"
target.models <- fread(file.path(problem.dir, "target_models.tsv"))

seg.file.dt <- data.table(seg.file=Sys.glob(file.path(problem.dir, "*_segments.bed")))

max.dt <- seg.file.dt[, {
  segs <- fread(seg.file)
  setnames(segs, c("chrom", "chromStart", "chromEnd", "status", "mean"))
  segs[, bases := chromEnd - chromStart]
  segs[, data.table(
    max.bases=max(bases),
    segments=.N
    )]
  }, by=seg.file]
setkey(max.dt, segments)
setkey(target.models, segments)
max.dt[target.models, .(fp, fn, segments, max.bases)]
