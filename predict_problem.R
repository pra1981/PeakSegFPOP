arg.vec <- c(#predict new data, target not yet computed.
  "test/H3K4me3_TDH_other/model.RData",
  "test/H3K4me3_TDH_other/samples/McGill0036/problems/24")
arg.vec <- c(#predicts inside, already computed.
  "test/H3K4me3_TDH_other/model.RData",
  "test/H3K4me3_TDH_other/samples/McGill0036/problems/22")
arg.vec <- c(#predicts outside target interval
  "test/H3K4me3_TDH_other/model.RData",
  "test/H3K4me3_TDH_other/samples/McGill0267/problems/7")
arg.vec <- c(#predicts infeasible.
  "test/H3K4me3_TDH_other/model.RData",
  "test/H3K4me3_TDH_other/samples/McGill0012/problems/24")

arg.vec <- commandArgs(trailingOnly=TRUE)

if(length(arg.vec) != 2){
  stop("usage: Rscript predict_problem.R model.RData problem_dir")
}
model.RData <- arg.vec[1]
problem.dir <- arg.vec[2]

library(data.table)
library(coseg)

problem.predict(problem.dir, model.RData)
