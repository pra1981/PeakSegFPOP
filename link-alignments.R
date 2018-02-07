for(experiment in c("H3K36me3", "H3K4me3")){
  for(types in c("immune", "other")){
    set.name <- paste0(experiment, "_TDH_", types)
    sample.dir.vec <- Sys.glob(file.path("labels", set.name, "samples", "*", "*"))
    sample.id.vec <- basename(sample.dir.vec)
    src.bed.vec <- file.path("~/bed", sample.id.vec, paste0(experiment, ".bed"))
    alignments.bed.vec <- file.path(sample.dir.vec, "alignments.bed")
    unlink(alignments.bed.vec)
    file.symlink(src.bed.vec, alignments.bed.vec)
  }
}
