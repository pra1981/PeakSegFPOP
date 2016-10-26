source("packages.R")

writeProblem <- function(data.list, problem.dir){
  file.list <- list(
    coverage=file.path(problem.dir, "coverage.bedGraph"),
    labels=file.path(problem.dir, "labels.bed"),
    problem=file.path(problem.dir, "problem.bed"))
  data.list$problem <- with(data.list$coverage, data.frame(
    chrom=chrom[1],
    chromStart=min(chromStart),
    chromEnd=max(chromEnd)))
  unlink(problem.dir, recursive=TRUE)
  dir.create(problem.dir, showWarnings=FALSE, recursive=TRUE)
  for(name in names(file.list)){
    df <- data.list[[name]]
    df$chromStart <- sprintf("%d", df$chromStart)
    df$chromEnd <- sprintf("%d", df$chromEnd)
    path <- file.list[[name]]
    write.table(
      df, path, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  }
}

download.to <- function
(u, f, writeFun=if(grepl("bigWig", f))writeBin else writeLines){
  if(!file.exists(f)){
    f.dir <- dirname(f)
    dir.create(f.dir, showWarnings=FALSE, recursive=TRUE)
    request <- GET(u)
    stop_for_status(request)
    writeFun(content(request), f)
  }
}
