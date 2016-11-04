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
label.txt <- "
chr10:33,061,897-33,162,814 noPeaks
chr10:33,168,111-33,196,481 peakStart bcell kidney
chr10:33,212,017-33,250,897 peakEnd bcell kidney
chr10:33,275,005-33,445,786 noPeaks
chr10:33,456,000-33,484,755 peakStart kidney
chr10:33,597,317-33,635,209 peakEnd kidney
chr10:33,662,034-33,974,942 noPeaks

chr10:35,182,820-35,261,001 noPeaks
chr10:35,261,418-35,314,654 peakStart bcell kidney
chr10:35,343,031-35,398,459 peakEnd bcell kidney

chr10:38,041,023-38,102,554 noPeaks
chr10:38,110,649-38,122,541 peakStart bcell kidney
chr10:38,138,274-38,191,035 peakEnd bcell kidney
chr10:38,191,040-38,221,846 noPeaks
chr10:38,226,848-38,242,263 peakStart bcell kidney
chr10:38,257,916-38,269,653 peakEnd bcell kidney
chr10:38,269,700-38,296,000 noPeaks
chr10:38,296,008-38,307,179 peakStart bcell kidney
chr10:38,379,045-38,391,967 peakStart bcell kidney
chr10:38,404,899-38,412,089 peakEnd bcell kidney
chr10:38,413,073-38,444,133 noPeaks

chr10:38,585,584-38,643,190 noPeaks
chr10:38,643,191-38,650,766 peakStart bcell kidney
chr10:38,731,066-38,750,574 peakEnd bcell kidney
chr10:38,750,960-38,790,663 noPeaks


"
set.dir <- file.path("test", "noinputlabels")
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
cmd <- paste("Rscript pipeline.R", set.dir)
status <- system(cmd)
test_that("pipeline script succeeds", {
  expect_equal(status, 0)
})
test_that("index.html is created", {
  expect_true(file.exists(file.path(set.dir, "index.html")))
})
