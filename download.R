bed <- read.table("http://hubs.hpc.mcgill.ca/~thocking/bed/index.txt", header=TRUE)
"http://hubs.hpc.mcgill.ca/~thocking/genomecov/H3K36me3/McGill0001.bedGraph"
pre <- "http://hubs.hpc.mcgill.ca/~thocking/genomecov/"

bed.file.vec <- paste(bed$file)
for(bed.file.i in seq_along(bed.file.vec)){
  bed.file <- bed.file.vec[[bed.file.i]]
  sample.id <- dirname(paste(bed.file))
  experiment <- sub(".bed$", "", basename(paste(bed.file)))
  bedGraph.file <- paste0(bed.file, "Graph")
  bedGraph.path <- file.path("data", bedGraph.file)
  if(!file.exists(bedGraph.path)){
    bed.url <- paste0(pre, experiment, "/", sample.id, ".bedGraph")
    cat(sprintf(
      "%4d / %4d %s -> %s\n",
      bed.file.i, length(bed.file.vec),
      bed.url, bedGraph.path))
    dir.create(dirname(bedGraph.path), showWarnings=FALSE, recursive=TRUE)
    download.file(bed.url, bedGraph.path)
  }
}

if(!file.exists("hg19.txt")){
  download.file("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/chromInfo.txt.gz", "hg19.txt.gz")
  system("gunzip hg19.txt.gz")
}

"grep chr10 H3K36me3.bed > H3K36me3_chr10.bed"
"genomeCoverageBed -bga -g ../../hg19.txt -i H3K36me3_chr10.bed > H3K36me3_chr10.bedGraph"
