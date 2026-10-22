# Run from the independent repository root.
stopifnot(file.exists("app.R"), file.exists("data/processed/visualization/manifest.json"),
          file.exists("fonts/NotoSerifCJKjp-Regular.otf"))
stage <- tempfile("kld-shinylive-app-")
dir.create(stage)
for (dir in c("scripts/visualization", "data/processed/visualization", "data/reference", "fonts")) {
  dir.create(file.path(stage, dir), recursive = TRUE, showWarnings = FALSE)
  files <- list.files(dir, full.names = TRUE)
  # Do not include source TSV duplicates or development files in the browser bundle.
  if (dir == "data/processed/visualization") files <- files[basename(files) == "manifest.json"]
  if (dir == "fonts") files <- files[basename(files) == "OFL.txt"]
  stopifnot(all(file.copy(files, file.path(stage, dir), overwrite = TRUE)))
}
file.copy("app.R", stage)
shinylive::export(stage, "_site", wasm_packages = TRUE,
                  template_params = list(title = "KLD — Interactive R plots"))
# Keep large data and fonts outside app.json; fetch only the selected table.
for (dir in c("data/processed/visualization", "fonts")) {
  dir.create(file.path("_site", dir), recursive = TRUE, showWarnings = FALSE)
  files <- list.files(dir, full.names = TRUE)
  files <- files[!grepl("\\.tsv$", files)]
  stopifnot(all(file.copy(files, file.path("_site", dir), overwrite = TRUE)))
}
file.create("_site/.nojekyll")
unlink(stage, recursive = TRUE)
