library(data.table)
readUCSC <- function(pre){
  pre.txt.gz <- paste0(pre, ".txt.gz")
  dl <- function(f){
    if(!file.exists(f)){
      u <- paste0(
        "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/",
        f)
      download.file(u, f)
    }
  }
  dl(pre.txt.gz)
  data.dt <- fread(paste("zcat", pre.txt.gz))
  pre.sql <- paste0(pre, ".sql")
  dl(pre.sql)
  sql.lines.vec <- readLines(pre.sql)
  names.sql.vec <- grep("^  `", sql.lines.vec, value=TRUE)
  names.dt <- fread(paste(names.sql.vec, collapse="\n"), sep="`", header=FALSE)
  setnames(data.dt, names.dt$V2)
  data.dt
}

refGene <- readUCSC("refGene")
##knownGene <- readUCSC("knownGene")
unique(refGene[, list(txStart)])
some.chr1 <- refGene[chrom=="chr1" & txEnd < 1e6]
some.chr1[, chromStart := txStart]
some.chr1[, chromEnd := txEnd]
cluster.dt <- data.table(PeakSegJoint::clusterPeaks(some.chr1))[, list(
  chromStart=min(chromStart),
  chromEnd=max(chromEnd),
  plus=sum(strand=="+"),#going right >>
  minus=sum(strand=="-")#going left <<<
  ), by=list(cluster)]
cluster.dt[, dist.to.next := c(chromStart[-1]-chromEnd[-.N], NA) ]

ggplot()+
  geom_segment(aes(
    txStart, seq_along(txStart),
    xend=txEnd, yend=seq_along(txEnd)),
               data=some.chr1)+
  geom_segment(aes(
    chromStart, 0, xend=chromEnd, yend=0), data=cluster.dt)
