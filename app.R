# Entry point shared by native Shiny and Shinylive.
source("scripts/visualization/browser_assets.R")
if (exists(".base_url", envir = .GlobalEnv, inherits = FALSE)) {
  .fetch_browser_asset("fonts/NotoSerifCJKjp-Regular.otf", "fonts/NotoSerifCJKjp-Regular.otf")
  dir.create("/home/web_user/fonts", recursive = TRUE, showWarnings = FALSE)
  file.copy("fonts/NotoSerifCJKjp-Regular.otf", "/home/web_user/fonts", overwrite = TRUE)
}
source("scripts/visualization/app.R", chdir = TRUE)$value
