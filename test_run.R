if(!require(testthat))install.packages("testthat")
if(!require(devtools))install.packages("devtools")
if(!require(coseg))devtools::install_github("tdhock/coseg")
if(!require(cosegData))devtools::install_github("tdhock/cosegData")
library(testthat)
getenv.or <- function(env.var, default){
  env.value <- Sys.getenv(env.var)
  if(env.value == ""){
    default
  }else{
    env.value
  }
}
file.name <- getenv.or("TEST_SUITE", "test_cases.R")
test_file(file.name, reporter=c("summary", "fail"))
