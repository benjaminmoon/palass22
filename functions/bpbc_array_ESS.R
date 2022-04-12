#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
# https://www.r-bloggers.com/2015/09/passing-arguments-to-an-r-script-from-command-lines/
if (length(args) != 1) {
  stop(
    "Exactly one argument must be supplied (full mcmc.log filename).\n",
    call. = FALSE
  )
}

library(coda)
library(data.table)

taxon <- args[1]
dirs <- "/work/glbcm/pg_pyrate/pyrate_mcmc_logs/"

for (i in seq_len(10)) {
  filename <-
    paste0(dirs, "Pg_", taxon, "_", i, taxon, "_seGibbsrj_mcmc.log")
  repname <-
    paste0(taxon, i)

  mcmc_log <-
    fread(filename)

  burn_in <- nrow(mcmc_log) * 0.25
  mcmc_log <-
    mcmc(
      mcmc_log,
      start = burn_in,
      thin  = 5000
    )

  png(
    paste(dirs, repname, "%d", "trace_plots.png", sep = "_"),
    height = 30,
    width  = 7,
    units  = "in",
    res    = 150
  )
  par(
    mfrow  = c(30, 1),
    mar    = c(2, 4, 2, 0),
    cex    = 0.5
  )
  traceplot(mcmc_log)
  dev.off()

  ess <- as.data.frame(effectiveSize(mcmc_log))
  colnames(ess) <- repname
  write.table(
    ess,
    file = paste0(dirs, repname, "_ess.tsv"),
    col.names = NA,
    row.names = TRUE
  )
}
