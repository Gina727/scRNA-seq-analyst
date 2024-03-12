library(plumber)
pr("api/plumber.R") %>%
  pr_run(port=8000)
