#!/bin/Rscript

suppressMessages(library(tidyverse))
options(scipen = 999)

parse_args <- function(args) {
  res_ <- list()
  i_ <- 1
  while (i_ <= length(args)) {
    if (startsWith(args[i_], "--")) {
      key_ <- sub("^--", "", args[i_])
      if (grepl("=", key_)) {
        parts <- strsplit(key_, "=")[[1]]
        res_[[parts[1]]] <- parts[2]
      } else if (i_ < length(args) && !startsWith(args[i_ + 1], "--")) {
        res_[[key_]] <- args[i_ + 1]
        i_ <- i_ + 1
      } else {
        res_[[key_]] <- TRUE  # boolean flag
      }
    }
    i_ <- i_ + 1
  }
  res_
}


n_digits_ <- function(x) {
  if (x == 0) return(1)
  floor(log10(abs(x))) + 1
}
