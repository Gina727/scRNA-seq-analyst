library("plumber")
pr("C:/Users/user/Documents/single_cell_rnaseq/app-api/app-api.R") %>%
# pr("api/app-api.R")
  pr_run(port = 8000)
