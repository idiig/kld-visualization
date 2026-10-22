# Shinylive binds the deployed webR URL in .GlobalEnv before loading the App.
# Use that URL so project Pages subpaths also resolve without a fixed hostname.
.fetch_browser_asset <- function(relative_path, dest) {
  if (file.exists(dest)) return(invisible(dest))
  if (!exists(".base_url", envir = .GlobalEnv, inherits = FALSE)) {
    stop("Missing local asset: ", dest)
  }
  base <- get(".base_url", envir = .GlobalEnv)
  stopifnot(grepl("/shinylive/webr/$", base))
  url <- paste0(sub("shinylive/webr/$", "", base), relative_path)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  ok <- FALSE
  on.exit(if (!ok) unlink(dest))
  download.file(url, dest, quiet = TRUE, mode = "wb")
  ok <- file.exists(dest) && file.size(dest) > 0
  if (!ok) stop("Unable to load asset: ", relative_path)
  invisible(dest)
}
