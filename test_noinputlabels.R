source("test_functions.R")

## Download bigWig files from github.
bigWig.part.vec <- c(
  "Input/MS010302",
  "bcell/MS010302",
  ## "Input/MS026601",
  ## "bcell/MS026601",
  ## "Input/MS002201",
  ## "kidney/MS002201",
  ## "Input/MS002202",
  "kidney/MS002202")
set.dir <- file.path("test", "input")
repos.url <- "https://raw.githubusercontent.com/tdhock/input-test-data/master/"
for(bigWig.part in bigWig.part.vec){
  bigWig.file <- file.path(set.dir, "samples", bigWig.part, "coverage.bigWig")
  bigWig.url <- paste0(repos.url, bigWig.part, ".bigwig")
  download.to(bigWig.url, bigWig.file)
}
labels.url <- paste0(repos.url, "kidney_bcell_labels.txt")
labels.file <- file.path(set.dir, "labels", "kidney_bcell_labels.txt")
download.to(labels.url, labels.file)
problems.bed <- file.path(set.dir, "problems.bed")
unlink(problems.bed)
system(paste("grep chr10 hg19_problems.bed | head >", problems.bed))

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
