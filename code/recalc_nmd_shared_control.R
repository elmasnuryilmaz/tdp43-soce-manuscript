#!/usr/bin/env Rscript
# Current NMD entrypoint: descriptive contrasts only. The previous four-contrast
# t/sign tests omitted shared-baseline covariance and must not be regenerated.
args <- commandArgs(trailingOnly = FALSE)
self <- sub("^--file=", "", args[grepl("^--file=", args)][1])
script <- file.path(dirname(normalizePath(self)), "audit_tables_v118.py")
python <- Sys.getenv("PYTHON", "python3")
status <- system2(python, shQuote(script))
quit(status=status)
