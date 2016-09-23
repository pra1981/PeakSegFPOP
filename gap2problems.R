## inputs: gap.bed chromInfo.txt

## outputs: problems.bed

arg.vec <- c("hg19_gap.bed", "hg19_chromInfo.txt", "hg19_problems.bed")
arg.vec <- commandArgs(trailingOnly=TRUE)
if(length(arg.vec) != 3){
  stop("usage: Rscript gap2problems.R gap.bed chromInfo.txt problems.bed")
}

library(data.table)

## data("hg19.gap", package="cosegData")
## write.table(
##   hg19.gap[,1:3], "hg19_gap.bed",
##   sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

gap.bed <- arg.vec[1]
chromInfo.txt <- arg.vec[2]
problems.bed <- arg.vec[3]

gap <- fread(gap.bed)
setnames(gap, c("chrom", "chromStart", "chromEnd"))
setkey(gap, chrom)

chromInfo <- fread(chromInfo.txt)
chromSizes <- chromInfo[, 1:2, with=FALSE]
setnames(chromSizes, c("chrom", "bases"))
setkey(chromSizes, chrom)

problems <- gap[, {
  bases <- chromSizes[chrom]$bases
  problemStart <- c(0, chromEnd)
  problemEnd <- c(chromStart, bases)
  data.table(problemStart, problemEnd)[problemStart < problemEnd,]
}, by=chrom]

write.table(problems, problems.bed, sep="\t", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
