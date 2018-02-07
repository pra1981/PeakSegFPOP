alignments.bam.vec <- Sys.glob("~/PeakSegFPOP/labels/*/samples/*/*/alignments.bam")
hg19_chromInfo.txt <- normalizePath("hg19_chromInfo.txt", mustWork=TRUE)

library(data.table)
for(i in seq_along(alignments.bam.vec)){
  alignments.bam <- alignments.bam.vec[[i]]
  cat(sprintf("%4d / %4d %s\n", i, length(alignments.bam.vec), alignments.bam))
  sdir <- dirname(alignments.bam)
  coverage.bigWig <- file.path(sdir, "coverage.bigWig")
  ##
  script.txt <- paste0("#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=24:00:00
#PBS -A bws-221-ae
#PBS -m ae
#PBS -M tdhock5@gmail.com
#PBS -V
#PBS -o ", coverage.bigWig, ".out
#PBS -e ", coverage.bigWig, ".err
#PBS -N cov", basename(dirname(coverage.bigWig)), "
bedtools genomecov -bg -g hg19_chromInfo.txt -ibam ", alignments.bam, " > ", coverage.bedGraph, "
bedGraphToBigWig ", coverage.bedGraph, " ", hg19_chromInfo.txt, " ", coverage.bigWig, "
rm ", coverage.bedGraph, "
")
  sh <- paste0(coverage.bigWig, ".sh")
  cat(script.txt, file=sh)
  if(!file.exists(coverage.bigWig)){
    system(paste("qsub", sh))
    coverage.bedGraph <- file.path(sdir, "coverage.bedGraph")
    bg.cmd <- paste(
      "bedtools genomecov -bg -g hg19_chromInfo.txt -ibam", 
      alignments.bam,
      ">",
      coverage.bedGraph)
    ##system(bg.cmd)
    cmd <- paste(
      "bedGraphToBigWig",
      coverage.bedGraph,
      "hg19_chromInfo.txt",
      coverage.bigWig)
    ##system(cmd)
    unlink(coverage.bedGraph)
  }
}
