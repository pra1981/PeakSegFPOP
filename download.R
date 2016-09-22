bed <- read.table("http://hubs.hpc.mcgill.ca/~thocking/bed/index.txt", header=TRUE)

"http://hubs.hpc.mcgill.ca/~thocking/genomecov/H3K36me3/McGill0001.bedGraph"

bed.file.vec <- bed$file[1]
for(bed.file.i in seq_along(bed.file.vec)){
  bed.file <- bed.file.vec[[bed.file.i]]
  bed.path <- file.path("data", bed.file)
  if(!file.exists(bed.path)){
    bed.url <- paste0(pre, bed.file)
    cat(sprintf(
      "%4d / %4d %s -> %s\n",
      bed.file.i, length(bed.file.vec),
      bed.url, bed.path))
    dir.create(dirname(bed.path), showWarnings=FALSE, recursive=TRUE)
    download.file(bed.url, bed.path)
  }
}

if(!file.exists("hg19.txt")){
  download.file("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/chromInfo.txt.gz", "hg19.txt.gz")
  system("gunzip hg19.txt.gz")
}

"grep chr10 H3K36me3.bed > H3K36me3_chr10.bed"
"genomeCoverageBed -bga -g ../../hg19.txt -i H3K36me3_chr10.bed > H3K36me3_chr10.bedGraph"
