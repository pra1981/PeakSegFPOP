library(coseg)
data("H3K4me3_XJ_immune_chunk1")
sample.id <- "McGill0106"
H3K4me3_XJ_immune_chunk1$count <- H3K4me3_XJ_immune_chunk1$coverage
by.sample <-
  split(H3K4me3_XJ_immune_chunk1, H3K4me3_XJ_immune_chunk1$sample.id)
one.sample <- by.sample[[sample.id]]
one.sample$chrom <- "chr1"
labels <- data.table(
  chrom="chr1",
  chromStart=20007000,
  chromEnd=20009000,
  annotation="peakStart")
sample.dir <- "labels/small/McGill0106"
dir.create(sample.dir, showWarnings=FALSE, recursive=TRUE)
coverage.bedGraph <- file.path(sample.dir, "coverage.bedGraph")
write.table(
  one.sample[, c("chrom", "chromStart", "chromEnd", "coverage")],
  coverage.bedGraph,
  quote=FALSE, sep="\t",
  row.names=FALSE, col.names=FALSE)
labels.bed <- file.path(sample.dir, "labels.bed")
write.table(
  labels, labels.bed,
  quote=FALSE, sep="\t",
  row.names=FALSE, col.names=FALSE)
