#!/usr/bin/env Rscript
# Regenerate every table and figure in the manuscript.
#
#   Rscript reproduce_paper.R                  # full run, auto-detects cores
#   Rscript reproduce_paper.R --quick          # 200 trials, fast check
#   Rscript reproduce_paper.R --nsim=2000 --cores=6 --out=results_v2

suppressPackageStartupMessages(library(adaptBOIN))

args  <- commandArgs(trailingOnly = TRUE)
nsim  <- 2000L
cores <- tryCatch(max(1L, parallel::detectCores(logical = TRUE) - 1L),
                  error = function(e) 1L)
out   <- "adaptboin_output"
figs  <- TRUE
study <- TRUE

for (a in args) {
  if      (grepl("^--nsim=",  a)) nsim  <- as.integer(sub("^--nsim=",  "", a))
  else if (grepl("^--cores=", a)) cores <- as.integer(sub("^--cores=", "", a))
  else if (grepl("^--out=",   a)) out   <- sub("^--out=", "", a)
  else if (a == "--quick")        nsim  <- 200L
  else if (a == "--nofig")        figs  <- FALSE
  else if (a == "--nostudy")      study <- FALSE
  else message("Ignoring unknown argument: ", a)
}

message(sprintf("nsim=%d  cores=%d  out=%s", nsim, cores, out))
reproduce_paper(out_dir = out, n_sim = nsim, n_cores = cores,
                do_figures = figs, do_studies = study)
