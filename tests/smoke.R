options(kld.plot.font_family = Sys.getenv("KLD_PLOT_FONT", "Noto Serif CJK JP"))
app <- source("app.R")$value
stopifnot(inherits(app, "shiny.appobj"), length(.plottable_tables(VIS_DIR)) == 8,
          identical(.infer_labels("entry-vs-gender-kokin.tsv")$value_label, "F"))
base <- file.path(VIS_DIR, "entry-vs-gender-kokin.tsv")
for (ext in c("png", "svg")) {
  out <- tempfile(fileext = paste0(".", ext))
  draw_point_vs_target(input = base, output = out, kld_variant = "norm", frequency_type = "term",
                       badge_selection = "frontier", frontier_n = 3,
                       value_label = "F", point_label = "entry", axis_label = "gender")
  stopifnot(file.size(out) > 1000)
  unlink(out)
}
# The former subprocess-only network panel must render with both file devices.
members <- data.frame(w1 = c("春", "花"), w2 = c("花", "月"),
                      assoc_w1_to_w2 = c(0.8, 0.2), assoc_w2_to_w1 = c(0.4, 0.7),
                      n_target = 3, n_not_target = 2)
for (ext in c("png", "svg")) {
  out <- tempfile(fileext = paste0(".", ext))
  .open_canvas(out, 450)
  result <- render_network_glyph_panel(members, "red", "F", "other", 0.01, 1)
  dev.off()
  stopifnot(file.size(out) > 1000, length(result$pending) == 0)
  unlink(out)
}
# Exercise the real Shiny action, table output and downloadable-render state.
shiny::testServer(server, {
  session$setInputs(input_file = "entry-vs-gender-kokin.tsv", input_mode = "cmdline",
                    cmdline_text = paste("input=data/processed/visualization/entry-vs-gender-kokin.tsv",
                                         "kld_variant=norm frequency_type=term badge_selection=frontier frontier_n=3"))
  session$setInputs(render = 1)
  stopifnot(file.exists(rendered()), !is.null(last_draw_args()), !is.null(annotated_table()))
  stopifnot(grepl("描画完了", status()))
})
cat("PASS: app, manifest, PNG/SVG, network PNG/SVG, Shiny render and table state\n")
