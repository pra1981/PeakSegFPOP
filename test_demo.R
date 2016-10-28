source("test_functions.R")

## Download bigWig files from github.
bigWig.part.vec <- c(
  "Input/MS010302",
  "bcell/MS010302",
  "Input/MS002202",
  "kidney/MS002202",
  "Input/MS026601",
  "bcell/MS026601",
  "Input/MS002201",
  "kidney/MS002201"
    )
label.txt <- "
chr10:33,061,897-33,162,814 noPeaks
chr10:33,456,000-33,484,755 peakStart kidney
chr10:33,597,317-33,635,209 peakEnd kidney
chr10:33,662,034-33,974,942 noPeaks

chr10:35,182,820-35,261,001 noPeaks
chr10:35,261,418-35,314,654 peakStart bcell kidney
chr10:35,343,031-35,398,459 peakEnd bcell kidney

chr10:38,041,023-38,102,554 noPeaks
chr10:38,296,008-38,307,179 peakStart bcell kidney
chr10:38,379,045-38,391,967 peakStart bcell kidney
chr10:38,404,899-38,412,089 peakEnd bcell kidney
chr10:38,413,073-38,444,133 noPeaks

chr10:38,585,584-38,643,190 noPeaks
chr10:38,643,191-38,650,766 peakStart bcell kidney
chr10:38,731,066-38,750,574 peakEnd bcell kidney
chr10:38,750,960-38,790,663 noPeaks

chr10:38,807,475-38,815,200 noPeaks
chr10:38,815,201-38,816,355 peakStart bcell kidney Input
chr10:38,818,377-38,819,342 peakEnd bcell kidney Input

chr10:39,098,319-39,111,384 noPeaks
chr10:39,125,134-39,125,550 peakStart bcell kidney Input
chr10:39,125,594-39,126,266 peakEnd bcell kidney Input
chr10:39,126,866-39,140,858 noPeaks
"
set.dir <- file.path("test", "demo")
repos.url <- "https://raw.githubusercontent.com/tdhock/input-test-data/master/"
for(bigWig.part in bigWig.part.vec){
  suffix <- ifelse(grepl("MS026601|MS002201", bigWig.part), "/", "_unlabeled/")
  bigWig.file <- file.path(set.dir, "samples", sub("/", suffix, bigWig.part), "coverage.bigWig")
  bigWig.url <- paste0(repos.url, bigWig.part, ".bigwig")
  download.to(bigWig.url, bigWig.file)
}
labels.file <- file.path(set.dir, "labels", "some_labels.txt")
dir.create(dirname(labels.file), showWarnings=FALSE, recursive=TRUE)
writeLines(label.txt, labels.file)
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

