bed <- read.table("http://hubs.hpc.mcgill.ca/~thocking/bed/index.txt", header=TRUE)
"http://hubs.hpc.mcgill.ca/~thocking/genomecov/H3K36me3/McGill0001.bedGraph"
pre <- "http://hubs.hpc.mcgill.ca/~thocking/genomecovgz/"

library(data.table)
load("labels.RData")
sample.ids <- labels[set.name=="H3K36me3_AM_immune", paste(unique(sample.id))]

bed.file.vec <- paste(subset(bed, !grepl("Input", file))$file)
bed.file.vec <- paste(subset(bed, grepl("H3K36me3", file))$file)
for(bed.file.i in seq_along(bed.file.vec)){
  bed.file <- bed.file.vec[[bed.file.i]]
  bedGraph.file <- paste0(bed.file, "Graph")
  bedGraph.path <- file.path("data", bedGraph.file)
  sample.id <- dirname(paste(bed.file))
  if(sample.id %in% sample.ids){
    experiment <- sub(".bed$", "", basename(paste(bed.file)))
    if(!file.exists(bedGraph.path)){
      bed.url <- paste0(pre, experiment, "/", sample.id, ".bedGraph.gz")
      gz.path <- paste0(bedGraph.path, ".gz")
      cat(sprintf(
        "%4d / %4d %s -> %s\n",
        bed.file.i, length(bed.file.vec),
        bed.url, gz.path))
      dir.create(dirname(gz.path), showWarnings=FALSE, recursive=TRUE)
      download.file(bed.url, gz.path)
      system(paste("gunzip", gz.path))
    }
  }else{
    cat(sprintf(
      "%4d / %4d rm %s\n",
      bed.file.i, length(bed.file.vec),
      bedGraph.path))
    unlink(bedGraph.path)
  }
}


if(!file.exists("hg19.txt")){
  download.file("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/chromInfo.txt.gz", "hg19.txt.gz")
  system("gunzip hg19.txt.gz")
}

"grep chr10 H3K36me3.bed > H3K36me3_chr10.bed"
"genomeCoverageBed -bga -g ../../hg19.txt -i H3K36me3_chr10.bed > H3K36me3_chr10.bedGraph"
