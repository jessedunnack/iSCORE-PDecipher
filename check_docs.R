#!/usr/bin/env Rscript
# Quick script to find undocumented functions

r_files <- list.files("/home/user/iSCORE-PDecipher/R", pattern = "\\.R$", full.names = TRUE)

undocumented <- list()

for (file in r_files) {
  lines <- readLines(file, warn = FALSE)

  for (i in seq_along(lines)) {
    # Check if line defines a function
    if (grepl("^[a-zA-Z._][a-zA-Z0-9._]* *<- *function", lines[i])) {
      # Check previous line for roxygen2 comment
      if (i > 1 && !grepl("^#'", lines[i-1])) {
        func_name <- sub(" *<-.*", "", lines[i])
        undocumented[[basename(file)]] <- c(undocumented[[basename(file)]],
                                           paste0("Line ", i, ": ", func_name))
      }
    }
  }
}

# Print results
for (file in names(undocumented)) {
  cat("\n===", file, "===\n")
  for (func in undocumented[[file]]) {
    cat("  ", func, "\n")
  }
}

cat("\n\nTotal undocumented functions:", sum(sapply(undocumented, length)), "\n")
