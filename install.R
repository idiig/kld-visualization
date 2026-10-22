packages <- c("shiny", "ggrepel", "data.table", "igraph", "jsonlite", "httpuv", "shinylive")
install.packages(packages, repos = "https://cloud.r-project.org",
                 Ncpus = max(1L, min(4L, parallel::detectCores())))
# install.packages() may warn and return success when a dependency failed.
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("R package installation failed: ", paste(missing, collapse = ", "))
