source("test_functions.R")
file.name <- getenv.or("TEST_SUITE", "test_cases.R")
test_file(file.name, reporter=c("summary", "fail"))
