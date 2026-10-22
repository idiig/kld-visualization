# Shiny replacement for dh-book-draft's scripts/analysis/demo/
# tune_pair_thresholds_server.py (a hand-rolled http.server tuning
# panel for plot_point_vs_target.R's own parameters) -- chosen over a
# plumber-based rewrite since draw_point_vs_target() is already a
# reusable sourced function (not CLI-only) and "reactive controls ->
# re-render plot" is exactly Shiny's own model (confirmed with the
# user, 2026-10-14).
#
# Deliberately a SMALLER scope than the original 1122-line tool, not a
# 1:1 port:
#   - lists whatever *.tsv ALREADY exists under data/processed/
#     visualization/ (produced by scripts/patches/visprepare.py +, for
#     pair tables, build_pair_enriched_table.R) -- does NOT auto-run
#     either of those when a file is missing (the original's own
#     auto-shell-out-to-demo/*.sh behavior). Regenerating a missing
#     slice is a manual step outside this app.
#   - covers draw_point_vs_target()'s own core/most-used parameters
#     (labels, badge_selection, aggregate_locators, edge_threshold/
#     dispersion_threshold, inset/zoom range -- added 2026-10-14 per
#     user request) -- not every CLI flag (network grid params,
#     annotate/highlight file uploads).
# Both are real, intentional reductions -- revisit if the user needs
# the dropped features back.
#
# Rendering is NOT automatic on every input change (unlike the
# original) -- draw_point_vs_target() can take real time on a 150k+row
# pair table, so re-rendering is gated behind an explicit "描画" button
# (actionButton + eventReactive), matching this app's own "don't
# surprise the user with a multi-second freeze on every keystroke"
# choice.
#
# Usage (works from any cwd, unlike every other Rscript in this
# project -- see below for why):
#   nix develop . -c Rscript -e 'shiny::runApp("scripts/visualization/app.R", port = 2718)'

library(shiny)

# Chromium blocks service-worker downloads when the anchor has a download
# attribute. Shiny supplies the filename in Content-Disposition already.
if (exists(".base_url", envir = .GlobalEnv, inherits = FALSE)) {
  downloadButton <- function(...) {
    tag <- shiny::downloadButton(...)
    tag$attribs$download <- NULL
    tag
  }
}

# shiny::runApp() on a single .R file setwd()s to that file's own
# directory before sourcing it (confirmed 2026-10-14 -- the earlier
# version's "scripts/visualization/plot_point_vs_target.R" relative
# path, written assuming cwd==repo root like every other Rscript
# invocation in this project, broke with "No such file or directory"
# under real runApp() because cwd is actually scripts/visualization/
# itself at this point, not the repo root). getwd() here IS reliably
# this file's own directory -- used to derive ROOT_DIR instead of
# assuming any particular caller cwd.
APP_DIR <- getwd()
ROOT_DIR <- normalizePath(file.path(APP_DIR, "..", ".."))

source(file.path(APP_DIR, "plot_point_vs_target.R"))
options(kld.plot.font_family = Sys.getenv("KLD_PLOT_FONT", "Noto Serif CJK JP"))
# Explicit imports also let Shinylive discover compiled package dependencies.
if (FALSE) {
  library(ggrepel)
  library(data.table)
  library(igraph)
  library(jsonlite)
}

VIS_DIR <- file.path(ROOT_DIR, "data", "processed", "visualization")
REF_DIR <- file.path(ROOT_DIR, "data", "reference")
.MANIFEST_PATH <- file.path(VIS_DIR, "manifest.json")
.DATA_MANIFEST <- if (file.exists(.MANIFEST_PATH)) jsonlite::read_json(.MANIFEST_PATH) else list()
.table_header <- function(path) {
  meta <- .DATA_MANIFEST[[basename(path)]]
  if (!is.null(meta)) return(unlist(meta$columns, use.names = FALSE))
  strsplit(readLines(path, n = 1, warn = FALSE), "\t")[[1]]
}

# renderImage() needs a real file at src -- before the first render
# (or after a failed one) there's no PNG yet, so output$preview falls
# back to this 1x1 transparent placeholder rather than an empty path,
# which renderImage would error on. Built once, module-level (not
# per-session) since it never changes.
.BLANK_PREVIEW_PNG <- local({
  path <- tempfile(fileext = ".png")
  png(path, width = 1, height = 1, bg = "transparent")
  dev.off()
  path
})

# Annotation-list files are picked from an existing directory (data/
# reference/), not uploaded -- same convention tune_pair_thresholds_
# server.py's own panel used ("this panel only picks which list file,
# it doesn't name the list itself"), chosen over a fileInput upload
# widget per the user's explicit choice (2026-10-14). "(none)" is
# always the first choice, mapping to NULL (the corresponding *_file
# argument stays unset).
.NONE <- "(none)"

.annotate_words_choices <- function(dir) {
  c(.NONE, list.files(dir, pattern = "\\.txt$", full.names = FALSE))
}

# highlight_file/highlight_pairs_file/annotate_pairs_file are all plain
# tsvs with a header row -- distinguished by their own column set, same
# three shapes resolve_annotation_inputs() itself expects:
#   highlight_file: term, type
#   highlight_pairs_file: w1, w2, type
#   annotate_pairs_file: w1, w2 (no type)
.tsv_choices_by_header <- function(dir, required_cols) {
  files <- list.files(dir, pattern = "\\.tsv$", full.names = FALSE)
  matching <- Filter(function(f) {
    header <- strsplit(readLines(file.path(dir, f), n = 1, warn = FALSE), "\t")[[1]]
    all(required_cols %in% header)
  }, files)
  c(.NONE, matching)
}

.resolve_ref_choice <- function(dir, choice) {
  if (identical(choice, .NONE) || is.null(choice) || choice == "") return(NULL)
  file.path(dir, choice)
}

# cmdlineモード(2026-10-17新增)用: ユーザーが貼り付けた生コマンド文字列
# ("再現コマンドをコピー"が吐くのと同じ形式、先頭の"Rscript ..."部分が
# あっても良い)からkey=valueトークンだけ拾って、plot_point_vs_target.R
# 側で定義済みの.kv_pairs()(source()で既に読み込まれている)に渡す --
# "="を含まないトークン(Rscriptそのものやスクリプトパス)は素通しで無視
# する。空文字/トークン無しならlength 0のcharacter vectorを返す。
.kv_from_cmdline <- function(text) {
  tokens <- strsplit(trimws(text), "\\s+")[[1]]
  tokens <- tokens[nzchar(tokens) & grepl("=", tokens, fixed = TRUE)]
  if (length(tokens) == 0) return(setNames(character(0), character(0)))
  .kv_pairs(tokens)
}

# input=/annotate_words_file=などのパス系の値はリポジトリルート相対で
# 書かれている前提("再現コマンドをコピー"自身がそう吐く形式) -- 絶対
# パスでなければROOT_DIRを前置する。
.resolve_cmdline_path <- function(path) {
  if (is.null(path) || !nzchar(path)) return(path)
  if (startsWith(path, "/")) path else file.path(ROOT_DIR, path)
}

# cmdlineモードで貼り付けられたkey=value(.kv_from_cmdlineの出力、パス系
# はまだリポジトリルート相対のまま、.resolve_cmdline_path前)を、対応する
# UI側のウィジェットにも反映する(2026-10-17新增、ユーザー要望"cmd里的
# 参数自动反映到ui") -- 「UIモードに切り替えた時に同じ設定から続けられる
# ように」が目的で、cmdlineモード自体の描画結果には影響しない(rendered()
# のcmdline分岐はinput$cmdline_textを毎回自分で読み直すので、このUI反映
# が無くても描画自体は正しく動く)。キー集合はALL_PARAM_KEYS + "input"
# (データセット選択だけ別input_idなので個別扱い)。コマンドに出てこな
# かったキーは、チェックボックス/フィルタ系は「オフ」、それ以外は各自の
# 元のデフォルト値に戻す -- 「貼り付けたコマンドの状態をそのまま映す」
# ことが目的なので、UI側に前から残っていた値が中途半端に混ざらないよう
# 毎回全部上書きする。
`%||%` <- function(a, b) if (is.null(a)) b else a

.sync_ui_from_kv <- function(session, kv) {
  .get <- function(key) if (key %in% names(kv)) kv[[key]] else NULL
  .bool <- function(key, default = FALSE) {
    v <- .get(key)
    if (is.null(v)) default else tolower(v) == "true"
  }
  .num <- function(key) {
    v <- .get(key)
    if (is.null(v)) NA_real_ else suppressWarnings(as.numeric(v))
  }
  .select_file <- function(id, key) {
    v <- .get(key)
    updateSelectInput(session, id, selected = if (is.null(v) || !nzchar(v)) .NONE else basename(v))
  }
  .sync_filter <- function(id) {
    v <- .num(id)
    if (!is.na(v)) {
      updateCheckboxInput(session, paste0("enable_", id), value = TRUE)
      updateNumericInput(session, id, value = v)
    } else {
      updateCheckboxInput(session, paste0("enable_", id), value = FALSE)
    }
  }

  if (!is.null(.get("input"))) {
    updateSelectInput(session, "input_file", selected = basename(.get("input")))
  }
  if (!is.null(.get("kld_variant"))) {
    updateSelectInput(session, "kld_variant", selected = .get("kld_variant"))
  }
  if (!is.null(.get("frequency_type"))) {
    updateSelectInput(session, "frequency_type", selected = .get("frequency_type"))
  }
  updateTextInput(session, "point_label", value = .get("point_label") %||% "point")
  updateTextInput(session, "axis_label", value = .get("axis_label") %||% "target")
  updateTextInput(session, "value_label", value = .get("value_label") %||% "value")
  if (!is.null(.get("badge_selection"))) {
    updateSelectInput(session, "badge_selection", selected = .get("badge_selection"))
  }
  updateCheckboxInput(session, "aggregate_locators", value = .bool("aggregate_locators", TRUE))

  for (id in c("edge_threshold", "dispersion_threshold", "node_dispersion_threshold",
               "frequency_threshold", "node_frequency_threshold")) {
    .sync_filter(id)
  }
  frontier_n <- .num("frontier_n")
  updateNumericInput(session, "frontier_n", value = if (is.na(frontier_n)) 15 else frontier_n)
  updateCheckboxInput(session, "annotate_nontarget_only", value = .bool("annotate_nontarget_only"))
  updateTextInput(session, "annotate_point_exclude", value = .get("annotate_point_exclude") %||% "")
  updateCheckboxInput(session, "highlight_frontier", value = .bool("highlight_frontier"))

  .select_file("annotate_words_file", "annotate_words_file")
  .select_file("highlight_file", "highlight_file")
  .select_file("highlight_pairs_file", "highlight_pairs_file")
  .select_file("annotate_pairs_file", "annotate_pairs_file")

  for (id in c("zoom_x1", "zoom_x2", "zoom_y1", "zoom_y2")) {
    updateNumericInput(session, id, value = .num(id))
  }
  updateCheckboxInput(session, "inset", value = .bool("inset"))
  if (!is.null(.get("inset_position"))) {
    updateSelectInput(session, "inset_position", selected = .get("inset_position"))
  }
  inset_pad_frac <- .num("inset_pad_frac")
  updateNumericInput(session, "inset_pad_frac", value = if (is.na(inset_pad_frac)) 0.15 else inset_pad_frac)
  updateCheckboxInput(session, "zoom_annotated_only", value = .bool("zoom_annotated_only"))
  updateCheckboxInput(session, "repel_annotated_only", value = .bool("repel_annotated_only"))
}

# point_label/axis_label/value_label aren't stored anywhere as English
# strings -- the calc yaml that produced a given visualization tsv only
# has structure.level1/level2 (Japanese field names, e.g. "和歌集"/
# "語"/"性別") and condition2.value (the actual target, e.g.
# "shinkokin"/"F"). This table translates the handful of field names
# this project's calc yamls actually use into the English labels
# already used throughout this session's own test renders -- not a
# fresh naming choice, just making that existing convention automatic.
# Display metadata exported independently of calculation configs.
.DATASET_LABELS <- list(
  "entry-vs-anthology-shinkokin.tsv" = list(point_label = "entry", axis_label = "anthology", value_label = "shinkokin"),
  "entry-vs-gender-first4.tsv" = list(point_label = "entry", axis_label = "gender", value_label = "F"),
  "entry-vs-gender-kokin.tsv" = list(point_label = "entry", axis_label = "gender", value_label = "F"),
  "pair-entry-entry-vs-anthology-shinkokin.tsv" = list(point_label = "entry", axis_label = "anthology", value_label = "shinkokin"),
  "pair-entry-entry-vs-gender-first4.tsv" = list(point_label = "entry", axis_label = "gender", value_label = "F"),
  "pair-entry-budatemajor-vs-gender-first4.tsv" = list(point_label = "budatemajor", axis_label = "gender", value_label = "F")
)

.infer_labels <- function(input_file, root_dir = ROOT_DIR) {
  .DATASET_LABELS[[input_file]]
}

.DEFAULT_INPUT_FILE <- "entry-vs-anthology-shinkokin.tsv"
.INITIAL_LABELS <- .infer_labels(.DEFAULT_INPUT_FILE)
if (is.null(.INITIAL_LABELS)) {
  .INITIAL_LABELS <- list(point_label = "point", axis_label = "target", value_label = "value")
}

# Only tsvs shaped like draw_point_vs_target()'s own input (a "point"
# header column) are offered -- association tables (w1/w2/assoc_w1_to_w2/
# assoc_w2_to_w1, no "point" column) are a different shape entirely and
# would just error if picked.
.plottable_tables <- function(dir) {
  files <- list.files(dir, pattern = "\\.tsv$", full.names = FALSE)
  files <- union(names(.DATA_MANIFEST), files)
  Filter(function(f) "point" %in% .table_header(file.path(dir, f)), files)
}

# An un-enriched pair table (has "w1" but not "freq1" -- never ran
# through build_pair_enriched_table.R) does NOT error in
# plot_point_vs_target.R -- verified directly (2026-10-14), it was
# assumed to crash in .aggregate_locator_groups (missing freq1/freq2/
# assoc_w1_to_w2/assoc_w2_to_w1 columns) but never actually reaches
# that code: .dedup_pair_rows (called first, for every badge_selection
# mode) does `rows$assoc_w1_to_w2` on a column that doesn't exist at
# all, getting NULL -> `!is.na(NULL)` is a zero-length logical ->
# `rows[logical(0), ]` silently selects ZERO rows regardless of
# NROW(rows). The real, user-visible effect is a plot with NO badges/
# labels at all, silently -- not a crash. Flagged here so the user
# doesn't mistake "no labels" for "nothing matched."
.is_unenriched_pair <- function(path) {
  header <- .table_header(path)
  "w1" %in% header && !("freq1" %in% header)
}

.pair_header_fields <- function(path) {
  .table_header(path)
}

# field1/field2 (when present at all -- a 2026-09-23 addition, see
# plot_point_vs_target.R's own .render_pair_network_svg docstring)
# name which of entry/concept/topic each side of the pair is. Every
# pair table this project actually has today predates that column (no
# entry+concept/entry+topic calc yaml exists yet -- concept2l2entry is
# fully unwired, see project memory) or simply never carries anything
# but entry+entry, so an absent field1/field2 is treated as entry+entry
# too, not "unknown."
.is_entry_entry_pair <- function(path, header) {
  if (!all(c("field1", "field2") %in% header)) return(TRUE)
  meta <- .DATA_MANIFEST[[basename(path)]]
  if (!is.null(meta)) return(identical(meta$field1, "entry") && identical(meta$field2, "entry"))
  row <- read.table(path, sep = "\t", header = TRUE, nrows = 1,
                     stringsAsFactors = FALSE, quote = "", fileEncoding = "UTF-8")
  identical(row$field1[1], "entry") && identical(row$field2[1], "entry")
}

# A small inline checkbox placed immediately beside its own option,
# toggling whether that option's key=value appears in the repro
# command printed at the bottom of the image (see .build_repro_cmd) --
# restores the old tool's own "per-row toggle, not one combined list"
# layout (tune_pair_thresholds_server.py's own injectCmdToggles, ~line
# 532), reversing this app's earlier single-checkboxGroupInput
# shortcut per the user's explicit request (2026-10-15). Lives inside
# the SAME conditionalPanel as its own option below -- a hidden
# option's toggle is hidden too, matching the old tool's own behavior
# (its own checkbox was injected onto the row itself, so row-hidden
# implied checkbox-hidden there too).
.cmd_toggle <- function(key) {
  tags$div(style = "margin: -6px 0 10px 18px; font-size: 0.85em; color: #666;",
           checkboxInput(paste0("cmd_incl_", key), "再現コマンドに含める", value = TRUE))
}

# フィルタ用——「有効にする」チェックボックスを入れるまで対応する
# numericInputは出ない(blank=無効という以前の規約から切り替え、2026-10-15
# ユーザー確認済み)。新しいfilterを1個増やすのはこの呼び出しを1行増やす
# だけ、という形にするためのヘルパー。PAIR_ONLYはedge_threshold/node_
# dispersion_threshold用——output.is_pair_inputでさらに1段くるむ。
.filter_control <- function(id, label, pair_only = FALSE) {
  body <- tagList(
    checkboxInput(paste0("enable_", id), paste0(label, "を有効にする"), value = FALSE),
    conditionalPanel(
      condition = paste0("input.enable_", id),
      numericInput(id, label, value = NA),
      .cmd_toggle(id)
    )
  )
  if (pair_only) conditionalPanel(condition = "output.is_pair_input", body) else body
}

# サーバー側の対応するthreshold解決——有効チェックが外れていれば(値が
# 残っていても)NULL、入っていてNAならNULL、それ以外はその値そのまま。
.resolve_threshold <- function(input, id) {
  if (!isTRUE(input[[paste0("enable_", id)]])) return(NULL)
  value <- input[[id]]
  if (is.na(value)) NULL else value
}

# Walked by .build_repro_cmd()'s own iteration below, and by server's
# own excluded-keys lookup (Filter(...,  input[[paste0("cmd_incl_", k)]])
# over this SAME list) -- one literal list, read from both places, so
# they can never drift apart (same reasoning tune_pair_thresholds_
# server.py's own CMD_FIELDS comment gives for deriving its list rather
# than hand-duplicating it, ~line 505-517).
ALL_PARAM_KEYS <- c(
  "kld_variant", "frequency_type", "point_label", "axis_label", "value_label", "badge_selection",
  "aggregate_locators", "edge_threshold", "dispersion_threshold",
  "node_dispersion_threshold", "frequency_threshold", "node_frequency_threshold", "frontier_n", "annotate_nontarget_only",
  "annotate_point_exclude", "highlight_frontier", "annotate_words_file", "highlight_file",
  "highlight_pairs_file", "annotate_pairs_file", "zoom_x1", "zoom_x2",
  "zoom_y1", "zoom_y2", "inset", "inset_position", "inset_pad_frac",
  "zoom_annotated_only", "repel_annotated_only"
)

# Pure: PARAMS (a named list, one entry per ALL_PARAM_KEYS -- already-
# NULL-resolved the same way rendered()'s own eventReactive resolves
# them, see server's own .resolve_render_params()) rendered as
# "Rscript scripts/visualization/plot_point_vs_target.R key=value ..."
# -- a NULL/FALSE-default/empty value is skipped the same way an
# unset CLI flag would be, and EXCLUDED_KEYS (the cmd_fields
# checkboxGroupInput's own unchecked items) is a second, independent
# reason to skip a key -- mirrors build_repro_cmd()'s own two
# conditions in the Python original (tune_pair_thresholds_server.py
# ~882-885): which flags render.R the IMAGE is never touched by
# EXCLUDED_KEYS, only which flags get PRINTED.
.build_repro_cmd <- function(input_path, params, excluded_keys) {
  kv <- Filter(Negate(is.null), lapply(ALL_PARAM_KEYS, function(key) {
    if (key %in% excluded_keys) return(NULL)
    value <- params[[key]]
    if (is.null(value) || identical(value, "") || identical(value, FALSE)) return(NULL)
    sprintf("%s=%s", key, paste(value, collapse = ","))
  }))
  paste(c(
    "Rscript scripts/visualization/plot_point_vs_target.R",
    sprintf("input=%s", input_path),
    unlist(kv)
  ), collapse = " ")
}

.repro_cmd_visual_width <- function(s) {
  sum(ifelse(utf8ToInt(enc2utf8(s)) > 0x7F, 2, 1))
}

# Pure: SVG_TEXT with one more text line appended at the bottom,
# containing CMD -- a direct R port of tune_pair_thresholds_server.py's
# own append_repro_cmd_line (~line 909): grows the root <svg>'s own
# height/viewBox by ROW_PT, then draws CMD as monospace text sized so
# its own visual width (CJK counted as width-2, same _visual_width
# heuristic) exactly fills the canvas -- not a fixed font-size, which
# would either overflow or look tiny depending on how long CMD is.
# Returns SVG_TEXT unchanged if the root tag doesn't match the expected
# shape (e.g. plot_point_vs_target.R's own svg() device output format
# changes) -- a failed post-process must never break the preview itself.
.append_repro_cmd_line <- function(svg_text, cmd, row_pt = 22) {
  if (!nzchar(cmd)) return(svg_text)
  pattern <- paste0(
    '(<svg[^>]*width="([0-9.]+)pt"[^>]*height=")([0-9.]+)',
    '(pt"[^>]*viewBox="[0-9.\\-]+ [0-9.\\-]+ [0-9.\\-]+ )([0-9.]+)(")'
  )
  m <- regmatches(svg_text, regexec(pattern, svg_text))[[1]]
  if (length(m) == 0) return(svg_text)
  width <- as.numeric(m[3])
  old_height <- as.numeric(m[4])
  new_height <- old_height + row_pt
  old_root <- m[1]
  new_root <- paste0(m[2], new_height, m[5], new_height, m[7])
  svg_text <- sub(old_root, new_root, svg_text, fixed = TRUE)
  font_size <- width / (max(1, .repro_cmd_visual_width(cmd)) * 0.6)
  text_el <- sprintf(
    '<text x="2" y="%f" font-family="monospace" font-size="%f">%s</text></svg>',
    old_height + row_pt * 0.7, font_size, htmltools::htmlEscape(cmd)
  )
  sub("</svg>\\s*$", text_el, svg_text)
}

ui <- fluidPage(
  # Shiny's client JS already toggles a "shiny-busy" class on <html>
  # for the whole time any computation (here: draw_point_vs_target(),
  # which can take real time on a 150k+row pair table -- see this
  # file's own top docstring) is outstanding on the server, from the
  # moment the request is sent until the response is flushed back --
  # the standard vanilla-Shiny way to show a busy indicator with no
  # extra package (shinycssloaders etc.) and no manual reactiveVal
  # flag-threading, since the class flips at the right moments
  # automatically. #generating-msg stays display:none the rest of the
  # time via this same rule.
  tags$head(tags$style(HTML(
    "#generating-msg { display: none; color: #c00; font-weight: bold; }
     html.shiny-busy #generating-msg { display: block; }"
  ))),
  titlePanel("point-vs-target 調整パネル"),
  sidebarLayout(
    sidebarPanel(
      # 2026-10-17新增: UI個別入力の代わりに、.build_repro_cmd()が吐く
      # ような"再現コマンド"をそのまま貼り付けて描画できるモード。cmdline
      # モードの時は下のUI一式(input_file選択からinset設定まで)を丸ごと
      # 隠す(ユーザー要望"cmd line 的情况下，不用展示其他的选项按钮")--
      # 貼り付けたコマンド文字列自体がinput=含めた全パラメータを持って
      # いるので、重複する個別入力を出す意味がない。
      radioButtons("input_mode", "入力方式",
                   choices = c("UI" = "ui", "コマンドライン" = "cmdline"), selected = "ui",
                   inline = TRUE),
      conditionalPanel(
        condition = "input.input_mode == 'ui'",
      selectInput("input_file", "入力ファイル (data/processed/visualization/)",
                  choices = .plottable_tables(VIS_DIR),
                  selected = .DEFAULT_INPUT_FILE),
      uiOutput("unenriched_warning"),
      selectInput("kld_variant", "kld_variant", choices = c("obs", "nofreq", "norm"),
                  selected = "norm"),
      .cmd_toggle("kld_variant"),
      selectInput("frequency_type", "frequency_type", choices = c("term", "doc"),
                  selected = "term"),
      .cmd_toggle("frequency_type"),
      # Initial values match .DEFAULT_INPUT_FILE's own .infer_labels()
      # result -- kept in sync with input_file's own "selected" above
      # by the observeEvent(input$input_file, ...) in server() below,
      # which re-infers on every dataset change, not just at startup.
      textInput("point_label", "point_label", .INITIAL_LABELS$point_label),
      .cmd_toggle("point_label"),
      textInput("axis_label", "axis_label", .INITIAL_LABELS$axis_label),
      .cmd_toggle("axis_label"),
      textInput("value_label", "value_label", .INITIAL_LABELS$value_label),
      .cmd_toggle("value_label"),
      selectInput("badge_selection", "badge_selection",
                  choices = c("words", "zoom_box", "frontier", "pairs"), selected = "words"),
      .cmd_toggle("badge_selection"),
      checkboxInput("aggregate_locators", "aggregate_locators", value = TRUE),
      .cmd_toggle("aggregate_locators"),
      tags$hr(),
      tags$b("フィルタ (filter)"),
      .filter_control("edge_threshold", "edge_threshold (pair only)", pair_only = TRUE),
      .filter_control("dispersion_threshold", "dispersion_threshold"),
      .filter_control("node_dispersion_threshold", "node_dispersion_threshold (pair片側)",
                      pair_only = TRUE),
      .filter_control("frequency_threshold", "frequency_threshold (total_freq)"),
      .filter_control("node_frequency_threshold", "node_frequency_threshold (pair片側freq1/freq2)",
                      pair_only = TRUE),
      conditionalPanel(
        condition = "input.badge_selection == 'frontier'",
        numericInput("frontier_n", "frontier_n (badge_selection=\"frontier\"時)", value = 15,
                     min = 1, step = 1),
        .cmd_toggle("frontier_n")
      ),
      checkboxInput("highlight_frontier",
                    "highlight_frontier (frequencyがfrontier中点より大きい点を染色、badge_selectionと独立)",
                    value = FALSE),
      .cmd_toggle("highlight_frontier"),
      conditionalPanel(
        condition = "input.badge_selection == 'words' && input.annotate_words_file != '(none)'",
        checkboxInput("annotate_nontarget_only",
                      "annotate_nontarget_only (n_target==0の候補のみbadge化)", value = FALSE),
        .cmd_toggle("annotate_nontarget_only"),
        textInput("annotate_point_exclude",
                 "annotate_point_exclude (除外するpoint文字列、カンマ区切り)", value = ""),
        .cmd_toggle("annotate_point_exclude")
      ),
      tags$hr(),
      tags$b("アノテーションリスト (data/reference/)"),
      conditionalPanel(
        condition = "input.badge_selection == 'words'",
        selectInput("annotate_words_file", "annotate_words_file",
                    choices = .annotate_words_choices(REF_DIR)),
        .cmd_toggle("annotate_words_file")
      ),
      conditionalPanel(
        condition = "!output.is_pair_input",
        selectInput("highlight_file", "highlight_file",
                    choices = .tsv_choices_by_header(REF_DIR, c("term", "type"))),
        .cmd_toggle("highlight_file")
      ),
      conditionalPanel(
        condition = "output.is_entry_entry_pair_input",
        selectInput("highlight_pairs_file", "highlight_pairs_file",
                    choices = .tsv_choices_by_header(REF_DIR, c("w1", "w2", "type"))),
        .cmd_toggle("highlight_pairs_file")
      ),
      conditionalPanel(
        condition = "input.badge_selection == 'pairs'",
        selectInput("annotate_pairs_file", "annotate_pairs_file (badge_selection=\"pairs\"時)",
                    choices = .tsv_choices_by_header(REF_DIR, c("w1", "w2"))),
        .cmd_toggle("annotate_pairs_file")
      ),
      tags$hr(),
      tags$b("範囲設定 (zoom_box)"),
      tags$p(style = "font-size: 0.85em; color: #666;",
             "badge_selection=\"zoom_box\"の取景条件、または拡大表示(inset)の",
             "拡大範囲として使われます。4つとも空欄ならデータのmin/max。"),
      conditionalPanel(
        condition = "input.badge_selection == 'zoom_box' || input.inset",
        fluidRow(
          column(6, numericInput("zoom_x1", "zoom_x1 (frequency 下限)", value = NA),
                 .cmd_toggle("zoom_x1")),
          column(6, numericInput("zoom_x2", "zoom_x2 (frequency 上限)", value = NA),
                 .cmd_toggle("zoom_x2"))
        ),
        fluidRow(
          column(6, numericInput("zoom_y1", "zoom_y1 (keyness 下限)", value = NA),
                 .cmd_toggle("zoom_y1")),
          column(6, numericInput("zoom_y2", "zoom_y2 (keyness 上限)", value = NA),
                 .cmd_toggle("zoom_y2"))
        )
      ),
      tags$hr(),
      checkboxInput("inset", "拡大表示 (inset)", value = FALSE),
      .cmd_toggle("inset"),
      conditionalPanel(
        condition = "input.inset",
        selectInput("inset_position", "inset_position",
                    choices = c("top", "bottom", "left", "right"), selected = "top"),
        .cmd_toggle("inset_position"),
        numericInput("inset_pad_frac", "inset_pad_frac", value = 0.15, step = 0.05),
        .cmd_toggle("inset_pad_frac"),
        checkboxInput("zoom_annotated_only", "zoom_annotated_only (拡大側は注釈された点のみ)",
                      value = FALSE),
        .cmd_toggle("zoom_annotated_only"),
        checkboxInput("repel_annotated_only", "repel_annotated_only (repel対象も注釈された点のみ)",
                      value = FALSE),
        .cmd_toggle("repel_annotated_only")
      )
      ),
      conditionalPanel(
        condition = "input.input_mode == 'cmdline'",
        tags$p(style = "font-size: 0.85em; color: #666;",
               "「再現コマンドをコピー」で出力されるのと同じ形式のkey=value",
               "を空白区切りで貼り付けてください(先頭の\"Rscript ...\"部分や",
               "output=はあっても無視されます、input=はdata/processed/",
               "visualization/配下のファイルをリポジトリルート相対で指定)。"),
        textAreaInput("cmdline_text", NULL, rows = 8, width = "100%",
                      placeholder = paste(
                        "input=data/processed/visualization/entry-vs-anthology-shinkokin.tsv",
                        "kld_variant=norm frequency_type=term badge_selection=words"
                      ))
      ),
      tags$hr(),
      actionButton("render", "描画", class = "btn-primary"),
      br(), br(),
      downloadButton("download_svg", "SVGを保存"),
      br(), br(),
      actionButton("copy_cmd", "再現コマンドをコピー")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("プロット",
          tags$div(id = "generating-msg", "画像生成中..."),
          # imageOutput()'s own default height is a fixed "400px" -- the
          # actual rendered PNG is usually taller than that (canvas_size
          # dependent), so the image overflowed its container and visually
          # overlapped whatever came after it (verbatimTextOutput("status"))
          # rather than pushing it down. height="auto" lets the container
          # grow to the image's own real size instead.
          imageOutput("preview", height = "auto"),
          verbatimTextOutput("status")
        ),
        # draw_point_vs_target()'s return value (see its own comment, 2026-10-18)
        # is the subset of the input tsv currently annotated on the plot --
        # same badge/frontier/zoom_box selection the image itself used, with
        # the real tsv columns (not the frequency/keyness axis-swap renames).
        tabPanel("表", tableOutput("annotated_table"))
      )
    )
  ),
  tags$script(HTML("
    Shiny.addCustomMessageHandler('copy_repro_cmd', function(cmd) {
      navigator.clipboard.writeText(cmd);
    });
  "))
)

server <- function(input, output, session) {
  output$unenriched_warning <- renderUI({
    path <- file.path(VIS_DIR, input$input_file)
    if (.is_unenriched_pair(path)) {
      tags$p(style = "color: darkorange;",
             "このデータにはラベル表示に必要な列が不足しています。",
             "エラーにはなりませんが、badge(ラベル)が1つも出ません",
             "必要な列を含むデータを選んでください。")
    }
  })

  output$is_pair_input <- reactive({
    "w1" %in% .pair_header_fields(file.path(VIS_DIR, input$input_file))
  })
  outputOptions(output, "is_pair_input", suspendWhenHidden = FALSE)

  output$is_entry_entry_pair_input <- reactive({
    path <- file.path(VIS_DIR, input$input_file)
    header <- .pair_header_fields(path)
    "w1" %in% header && .is_entry_entry_pair(path, header)
  })
  outputOptions(output, "is_entry_entry_pair_input", suspendWhenHidden = FALSE)

  # Re-infer point_label/axis_label/value_label from the newly-picked
  # dataset's own calc yaml every time input_file changes (not just at
  # startup, see .INITIAL_LABELS above for the startup copy) --
  # unconditionally overwrites whatever is currently in the three text
  # inputs, including a value the user typed by hand, since a label
  # left over from the PREVIOUS dataset is normally just wrong for the
  # new one. ignoreInit=TRUE -- the startup value already came from
  # .INITIAL_LABELS, no need to redo the same lookup on first load.
  observeEvent(input$input_file, {
    labels <- .infer_labels(input$input_file)
    if (is.null(labels)) return(invisible(NULL))
    updateTextInput(session, "point_label", value = labels$point_label)
    updateTextInput(session, "axis_label", value = labels$axis_label)
    updateTextInput(session, "value_label", value = labels$value_label)
  }, ignoreInit = TRUE)

  # Holds the exact draw_point_vs_target() kwargs behind the last
  # successful render -- download_svg's own content() reuses these,
  # swapping only OUTPUT, to lazily regenerate a real vector file on
  # demand (see rendered()'s own comment for why the live preview no
  # longer renders svg at all).
  last_draw_args <- reactiveVal(NULL)
  status <- reactiveVal("")
  repro_cmd <- reactiveVal("")
  annotated_table <- reactiveVal(NULL)

  # cmdlineモードで貼り付けた内容を、同じ「描画」クリックのタイミングで
  # UI側のウィジェットにも反映する(.sync_ui_from_kv自身の docstring参照)
  # -- rendered()とは別のobserveEventだが、同じinput$renderイベントに
  # 両方ぶら下げているだけなので、描画自体の成否には影響しない。
  observeEvent(input$render, {
    if (!identical(input$input_mode, "cmdline")) return(invisible(NULL))
    .sync_ui_from_kv(session, .kv_from_cmdline(input$cmdline_text))
  })

  rendered <- eventReactive(input$render, {
    # PNG, not svg -- a pair dataset's scatter can be ~300k points,
    # which as svg() vector elements is both slow to write and (once
    # embedded as inline HTML) slow for a browser to render; a raster
    # device's cost scales with pixel count, not point count (see
    # .open_canvas's own comment in plot_point_vs_target.R, 2026-10-17).
    # The real vector .svg is only ever generated lazily, inside
    # download_svg's own content() below, when a user actually wants
    # one to put in the book.
    out <- tempfile(fileext = ".png")

    # cmdlineモード: UI側の個別入力は一切見ず、貼り付けられた生コマンド
    # だけからdraw_point_vs_target()の引数一式を組み立てる -- CLI入口
    # (plot_point_vs_target.Rの.draw_point_vs_target_args())と同じ
    # マッピングを再利用するので、どのパラメータがどう解釈されるかは
    # 普通にRscriptで実行した場合と完全に一致する。
    if (identical(input$input_mode, "cmdline")) {
      kv <- .kv_from_cmdline(input$cmdline_text)
      if (length(kv) == 0 || !("input" %in% names(kv))) {
        status("エラー: input=... を含むコマンドを貼り付けてください")
        return(NULL)
      }
      kv[["input"]] <- .resolve_cmdline_path(kv[["input"]])
      for (path_key in c("annotate_words_file", "highlight_file",
                          "highlight_pairs_file", "annotate_pairs_file")) {
        if (path_key %in% names(kv)) kv[[path_key]] <- .resolve_cmdline_path(kv[[path_key]])
      }
      draw_args <- .draw_point_vs_target_args(kv, output = out)
      ok <- tryCatch({
        annotated_table(do.call(draw_point_vs_target, draw_args))
        TRUE
      }, error = function(e) {
        status(paste("エラー:", conditionMessage(e)))
        FALSE
      })
      if (!ok || !file.exists(out)) return(NULL)
      status(paste("描画完了 (cmdline):", basename(kv[["input"]])))
      last_draw_args(draw_args)
      repro_cmd(paste(
        "Rscript scripts/visualization/plot_point_vs_target.R",
        paste(names(kv), kv, sep = "=", collapse = " ")
      ))
      return(out)
    }

    input_path <- file.path(VIS_DIR, input$input_file)
    edge_threshold <- .resolve_threshold(input, "edge_threshold")
    dispersion_threshold <- .resolve_threshold(input, "dispersion_threshold")
    node_dispersion_threshold <- .resolve_threshold(input, "node_dispersion_threshold")
    frequency_threshold <- .resolve_threshold(input, "frequency_threshold")
    node_frequency_threshold <- .resolve_threshold(input, "node_frequency_threshold")
    zoom_x1 <- if (is.na(input$zoom_x1)) NULL else input$zoom_x1
    zoom_x2 <- if (is.na(input$zoom_x2)) NULL else input$zoom_x2
    zoom_y1 <- if (is.na(input$zoom_y1)) NULL else input$zoom_y1
    zoom_y2 <- if (is.na(input$zoom_y2)) NULL else input$zoom_y2
    exclude_raw <- trimws(input$annotate_point_exclude)
    annotate_point_exclude <- if (exclude_raw == "") NULL else strsplit(exclude_raw, ",", fixed = TRUE)[[1]]
    annotate_words_file <- .resolve_ref_choice(REF_DIR, input$annotate_words_file)
    highlight_file <- .resolve_ref_choice(REF_DIR, input$highlight_file)
    highlight_pairs_file <- .resolve_ref_choice(REF_DIR, input$highlight_pairs_file)
    annotate_pairs_file <- .resolve_ref_choice(REF_DIR, input$annotate_pairs_file)
    draw_args <- list(
      input = input_path, output = out,
      kld_variant = input$kld_variant,
      frequency_type = input$frequency_type,
      point_label = input$point_label, axis_label = input$axis_label,
      value_label = input$value_label,
      badge_selection = input$badge_selection,
      aggregate_locators = input$aggregate_locators,
      edge_threshold = edge_threshold, dispersion_threshold = dispersion_threshold,
      node_dispersion_threshold = node_dispersion_threshold,
      frequency_threshold = frequency_threshold,
      node_frequency_threshold = node_frequency_threshold,
      frontier_n = input$frontier_n,
      annotate_nontarget_only = input$annotate_nontarget_only,
      annotate_point_exclude = annotate_point_exclude,
      highlight_frontier = input$highlight_frontier,
      annotate_words_file = annotate_words_file,
      highlight_file = highlight_file,
      highlight_pairs_file = highlight_pairs_file,
      annotate_pairs_file = annotate_pairs_file,
      zoom_x1 = zoom_x1, zoom_x2 = zoom_x2, zoom_y1 = zoom_y1, zoom_y2 = zoom_y2,
      inset = input$inset, inset_position = input$inset_position,
      inset_pad_frac = input$inset_pad_frac,
      zoom_annotated_only = input$zoom_annotated_only,
      repel_annotated_only = input$repel_annotated_only
    )
    # Success/failure is read off whether it threw AND whether OUT
    # actually got written, not off the return value -- but the return
    # value itself (the currently-annotated points, see draw_point_vs_
    # target()'s own 2026-10-18 comment) is still captured for the table
    # tab.
    ok <- tryCatch({
      annotated_table(do.call(draw_point_vs_target, draw_args))
      TRUE
    }, error = function(e) {
      status(paste("エラー:", conditionMessage(e)))
      FALSE
    })
    if (!ok || !file.exists(out)) return(NULL)
    status(paste("描画完了:", input$input_file))
    last_draw_args(draw_args)

    # Repro command: same resolved values as the real draw_point_vs_
    # target() call above, EXCEPT the four file fields print a repo-
    # root-relative path (data/reference/<file>) instead of the
    # absolute REF_DIR path this app itself resolved them to -- a
    # pasted command runs from repo root (same convention every other
    # Rscript invocation in this project already follows), not from
    # wherever this app's own REF_DIR/VIS_DIR constants happen to
    # point on this machine.
    .rel_ref <- function(choice) if (identical(choice, .NONE)) NULL else file.path("data", "reference", choice)
    cmd_params <- list(
      kld_variant = input$kld_variant,
      frequency_type = input$frequency_type,
      point_label = input$point_label, axis_label = input$axis_label,
      value_label = input$value_label, badge_selection = input$badge_selection,
      aggregate_locators = input$aggregate_locators,
      edge_threshold = edge_threshold, dispersion_threshold = dispersion_threshold,
      node_dispersion_threshold = node_dispersion_threshold,
      frequency_threshold = frequency_threshold,
      node_frequency_threshold = node_frequency_threshold,
      frontier_n = input$frontier_n,
      annotate_nontarget_only = input$annotate_nontarget_only,
      annotate_point_exclude = annotate_point_exclude,
      highlight_frontier = input$highlight_frontier,
      annotate_words_file = .rel_ref(input$annotate_words_file),
      highlight_file = .rel_ref(input$highlight_file),
      highlight_pairs_file = .rel_ref(input$highlight_pairs_file),
      annotate_pairs_file = .rel_ref(input$annotate_pairs_file),
      zoom_x1 = zoom_x1, zoom_x2 = zoom_x2, zoom_y1 = zoom_y1, zoom_y2 = zoom_y2,
      inset = input$inset, inset_position = input$inset_position,
      inset_pad_frac = input$inset_pad_frac,
      zoom_annotated_only = input$zoom_annotated_only,
      repel_annotated_only = input$repel_annotated_only
    )
    excluded_keys <- Filter(
      function(k) isFALSE(input[[paste0("cmd_incl_", k)]]), ALL_PARAM_KEYS
    )
    cmd <- .build_repro_cmd(
      file.path("data", "processed", "visualization", input$input_file),
      cmd_params, excluded_keys
    )
    repro_cmd(cmd)
    out
  })

  output$preview <- renderImage({
    path <- rendered()
    if (is.null(path)) return(list(src = .BLANK_PREVIEW_PNG, contentType = "image/png"))
    list(src = path, contentType = "image/png", style = "max-width: 100%; height: auto;")
  }, deleteFile = FALSE)

  output$annotated_table <- renderTable({
    tbl <- annotated_table()
    if (is.null(tbl) || nrow(tbl) == 0) return(data.frame())
    tbl
  })

  output$status <- renderText(status())

  output$download_svg <- downloadHandler(
    filename = function() paste0(tools::file_path_sans_ext(input$input_file), ".svg"),
    content = function(file) {
      # The live preview only ever renders png (see rendered()'s own
      # comment) -- the real vector svg is generated here, lazily, on
      # the SAME params the last successful preview used, only when a
      # user actually wants one to put in the book.
      draw_args <- last_draw_args()
      if (is.null(draw_args)) stop("先に「描画」を押してください")
      svg_out <- tempfile(fileext = ".svg")
      draw_args$output <- svg_out
      do.call(draw_point_vs_target, draw_args)
      writeLines(
        .append_repro_cmd_line(paste(readLines(svg_out, warn = FALSE), collapse = "\n"), repro_cmd()),
        svg_out
      )
      file.copy(svg_out, file)
    }
  )

  observeEvent(input$copy_cmd, {
    session$sendCustomMessage("copy_repro_cmd", repro_cmd())
  })
}

shinyApp(ui, server)
