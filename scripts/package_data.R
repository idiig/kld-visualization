# Serialization only: no statistical calculation or filtering.
# Usage: Rscript scripts/package_data.R /path/to/precomputed/visualization
source("scripts/visualization/plot_point_vs_target.R")
args <- commandArgs(TRUE)
src <- if (length(args)) args[[1]] else "data/processed/visualization"
dest <- "data/processed/visualization"
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
manifest <- list()
for (path in sort(list.files(src, pattern = "\\.tsv$", full.names = TRUE))) {
  header <- strsplit(readLines(path, n = 1, warn = FALSE), "\t")[[1]]
  if (!identical(header[[1]], "point")) next
  d <- load_and_merge_input(path)
  name <- basename(path)
  asset <- sub("\\.tsv$", ".rds", name)
  saveRDS(d, file.path(dest, asset), compress = "gzip", version = 3)
  stopifnot(identical(d, readRDS(file.path(dest, asset))))
  manifest[[name]] <- list(asset = asset, rows = nrow(d), columns = names(d),
                         field1 = if ("field1" %in% names(d)) d$field1[[1]] else NULL,
                         field2 = if ("field2" %in% names(d)) d$field2[[1]] else NULL,
                         source_md5 = unname(tools::md5sum(path)))
  message(name, ": ", nrow(d), " rows, ", round(file.size(file.path(dest, asset))/1024^2, 2), " MiB")
}
stopifnot(length(manifest) > 0)
jsonlite::write_json(manifest, file.path(dest, "manifest.json"), pretty = TRUE, auto_unbox = TRUE)
