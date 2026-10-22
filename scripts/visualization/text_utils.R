stem_tokens_ge2 <- function(stem_col) {
  # 单字助词/助动词词干在entry字段里到处撞车，只保留长度>=2的词根
  stems <- trimws(unlist(strsplit(stem_col, "\\|")))
  unique(stems[nchar(stems) >= 2])
}
