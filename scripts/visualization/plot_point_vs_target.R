# 可复用的"点 vs 固定值轴"散点图，风格照搬scripts/kld/plot_pair_frequency_keyness.R/
# plot_dispersion_extreme50_scatter.R(进而是scripts/hapax_male.R)已经定下来的一套
# 约定，不是另起一套：
#   - 不画box()外框，只用plot.new()+plot.window()+axis(1)/axis(2)画左/下轴线
#   - 正方形画布(2026-09-13从480pt高的黄金比例长方形改过来，边长后来又从360
#     改到450)，canvas_scale按画布边长/300pt换算mar/mgp/cex
#   - 坐标轴标题用bquote()做KLD_signed数学排版，数值量级小时刻度除以一个
#     比例显示、比例写进标题括号里(不是让R自己出"5.0e-06"这种科学计数法)
#   - 两轴都画abline(h=neutral,v=neutral)虚线(理论中性点：raw模式是0，
#     nofreq模式是0.5)，另加一条中位数实线(纯经验对半分，不跨图可比，不跟
#     虚线混用同一种线型)
#   - 轴两端标value_label/not_value_label代替看正负号，放在所有点画完之后
#     (最后画、盖在最上层)，x轴右端和y轴顶端各标一次value_label，x轴左端/
#     y轴底端合并成一个not_value_label(不然会在左下角重叠两次)
#
# x/y轴含义: x=frequency(TARGET measured against POINT，"TARGET→POINT"方向，
# corpus→word)，y=keyness(POINT measured against TARGET，"POINT→TARGET"方向，
# word→corpus)——两个轴都是同一个relationship的signed KLD读数，不是"一个KLD、
# 一个粗频度"。
#
# 输入表(见extract_point_vs_target.py)：列 point, keyness, frequency, total_freq,
# n_target, n_not_target, dispersion。
#
# 用法(key=value风格，跟kld_nofreq自己的CLI约定一致):
#   Rscript scripts/analysis/plot_point_vs_target.R \
#     input=data/processed/point-vs-target/entry-vs-gender-combined.tsv \
#     output=artifacts/figs/entry-vs-gender-combined.svg \
#     point_label=語 axis_label=ジェンダ value_label=女
#
# 也可以被run_meta_config_plots.R直接source()、调用draw_point_vs_target()
# 函数复用，不必每次都开一个新的Rscript进程。
#
# 2026-09-13模块化重构：原来一个约420行的函数靠嵌套closure分工，现在拆成
# 下面这些顶层具名函数，each职责单一，draw_point_vs_target()本身只负责按
# 顺序调度，不包含具体绘图细节。同时把两个写死语义的参数(label_words=badge、
# gradient_words=染色)通用化成两份独立命名的名单：annotate_words(单列，
# 画badge)、highlight(term+type两列，按type分类别染色)，支持任意数量的
# type，每个type自动分配一条渐变色带、进图例；badge排版加了第二种可选
# 算法(最大空白矩形，新默认)，跟原有的力导向repel(重新引入
# ggrepel)并列，用badge_layout参数切换。

library(ggrepel)

# ---- 坐标轴/画布 ----

# nofreq把signed KLD线性映射到[0,1](kld_lower->0, kld_upper->1，见
# kld.py:kld_and_nofreq的docstring)，不是以0为中心的signed范围——0本身
# 是"全部落在other一侧"的极端值，不是"无关联"的中性点，0.5(量表中点)才是
# 该画参考虚线的位置。raw(signed)模式保持"以0为中性点、按数据范围留白"的
# 画法。nofreq轴范围基本上钉在[0,1]，但留一点点边距(3%)——median(keyness)/
# median(frequency)经常恰好等于0或1(样本本身就这么分布，2026-09-14在kokin
# 上实测到median就是0)，完全不留白的话median那条辅助线会正好画在坐标轴
# 边框上，被边框本身盖住看起来像没画出来；留一点点边距能让贴边的辅助线
# 可见，同时因为3%很小，看起来仍然跟"钉死在[0,1]"没有明显区别。
.compute_axis_ranges <- function(d, nofreq) {
  if (nofreq) {
    list(xr = c(0, 1), yr = c(0, 1), xpad = 0.03, ypad = 0.03, neutral = 0.5)
  } else {
    xr <- range(d$frequency)
    yr <- range(d$keyness)
    # 留白比例从8%调到15%——跟nofreq模式现在的视觉留白比对下偏小
    # (2026-09-14确认"两者的视觉留白需要近似，目前非nofreq的留白过小")
    list(xr = xr, yr = yr, xpad = diff(xr) * 0.15, ypad = diff(yr) * 0.15, neutral = 0)
  }
}

# 参考scripts/hapax_pairs_male.R的"//"断轴标记做法：nofreq理论上限是1，但
# 实际数据经常远小于1(比如frequency方向可能只到0.03)——钉死在[0,1]会让
# 全部点挤在坐标轴一小段里。数据最大值离理论上限的差距超过坐标轴总范围
# 的break_frac(默认50%)时，把[0, axis_max]拆成两段：[0, main_end](数据
# 实际所在范围，留15%余量)线性映射到显示范围的前main_frac(默认82%)，
# [main_end, axis_max]压缩映射到剩下的部分，两段之间画"//"断轴标记
# (2026-09-14确认"两个值相差这个坐标轴大小超过50%的时候做轴的省略")。
# 返回NULL表示这根轴不需要断轴(数据本来就占满了理论范围，正常线性画)。
.axis_break_transform <- function(data_max, axis_max = 1, break_frac = 0.5, main_frac = 0.82) {
  if (!is.finite(data_max) || data_max <= 0) return(NULL)
  gap <- axis_max - data_max
  if (gap / axis_max < break_frac) return(NULL)
  main_end <- min(data_max * 1.15, axis_max * 0.999)
  fwd <- function(x) {
    ifelse(x <= main_end,
           x / main_end * main_frac,
           main_frac + (x - main_end) / (axis_max - main_end) * (1 - main_frac))
  }
  list(fwd = fwd, main_end = main_end, main_frac = main_frac, axis_max = axis_max)
}

# nofreq模式下检测x(frequency)/y(keyness)两根轴是否需要断轴，需要的话把
# d对应列原地替换成变换后的坐标——下游(点的位置、badge排版、连线、
# abline的中位数线)全部继续读d$frequency/d$keyness，不需要额外改动，
# 因为它们看到的就已经是"压缩过的"坐标；只有画刻度/断轴标记/neutral参考
# 线需要知道有没有断轴、断在哪，所以额外返回x_break/y_break。raw(signed)
# 模式没有固定的理论上限，不做断轴处理。
.apply_nofreq_breaks <- function(d, nofreq) {
  if (!nofreq) return(list(d = d, x_break = NULL, y_break = NULL))
  x_break <- .axis_break_transform(max(d$frequency))
  y_break <- .axis_break_transform(max(d$keyness))
  if (!is.null(x_break)) d$frequency <- x_break$fwd(d$frequency)
  if (!is.null(y_break)) d$keyness <- y_break$fwd(d$keyness)
  list(d = d, x_break = x_break, y_break = y_break)
}

# 在坐标轴基线上画一个"//"断轴标记，pos_frac是断点在[0,1]显示范围里的
# 相对位置，side=1画在x轴(底边)上，side=2画在y轴(左边)上——跟
# scripts/hapax_pairs_male.R的segments()画法同一个视觉(两条穿过轴线的
# 短斜线)，只是那边是按离散bar位置手算，这里按usr直接算，更通用。
.draw_axis_break_mark <- function(pos_frac, side, usr) {
  if (side == 1) {
    xm <- usr[1] + pos_frac * (usr[2] - usr[1])
    y0 <- usr[3]
    dy <- (usr[4] - usr[3]) * 0.018
    dx <- (usr[2] - usr[1]) * 0.008
    segments(xm - dx, y0 - dy, xm + dx, y0 + dy, xpd = NA, lwd = 1.5)
    segments(xm + dx, y0 - dy, xm + 3 * dx, y0 + dy, xpd = NA, lwd = 1.5)
  } else {
    ym <- usr[3] + pos_frac * (usr[4] - usr[3])
    x0 <- usr[1]
    dx <- (usr[2] - usr[1]) * 0.018
    dy <- (usr[4] - usr[3]) * 0.008
    segments(x0 - dx, ym - dy, x0 + dx, ym + dy, xpd = NA, lwd = 1.5)
    segments(x0 - dx, ym + dy, x0 + dx, ym + 3 * dy, xpd = NA, lwd = 1.5)
  }
}

# 断轴的刻度：主段(0到main_end)用pretty()正常取整数刻度，额外补一个
# axis_max(理论上限，比如nofreq的1)本身，每个刻度值都用fwd()换算成显示
# 位置，标签仍然是原始(未压缩)数值——这样读者看到的刻度值是真实的，只是
# 刻度间距不再均匀(主段疏、压缩段的最后一个点单独挤在末尾)。
.axis_ticks_with_break <- function(brk) {
  ticks <- pretty(c(0, brk$main_end))
  ticks <- ticks[ticks <= brk$main_end]
  ticks <- unique(c(ticks, brk$axis_max))
  list(at = brk$fwd(ticks), labels = format(ticks, trim = TRUE))
}

# width_multiplier/height_multiplier给inset模式用：画布整体加宽或加高一倍
# (另一边不变)，两半各自正好是一个canvas_size正方形(用par(fig=)切出来，
# 见.render_with_inset；2026-09-21改成上下叠放前用的是width_multiplier=2
# 左右分栏，现在默认改用height_multiplier=2上下分栏，width_multiplier保留
# 但目前没有调用方在用)。mar/mgp/cex都只由canvas_size决定，不受
# width_multiplier/height_multiplier影响——保证每一半看起来都跟一张独立的
# canvas_size图完全一样，不会因为共享一个更宽/更高的画布而字体/边距比例
# 跑掉。
.open_canvas <- function(output, canvas_size, width_multiplier = 1, height_multiplier = 1) {
  dir.create(dirname(output), showWarnings = FALSE, recursive = TRUE)
  canvas_scale <- canvas_size / 300
  panel_h_in <- canvas_size * height_multiplier / 72
  width_in <- canvas_size * width_multiplier / 72
  # PNG output -- a raster device's file size/render cost scales with
  # PIXEL count, not point count, unlike svg()'s one vector element
  # per point. A pair dataset's ~300k-point scatter makes svg() both
  # slow to write and (once embedded inline, see the Shiny app's own
  # output$preview) slow for a browser to parse -- confirmed 2026-10-17,
  # a 300k-point render was a 100MB+ SVG. This project's own LaTeX
  # pipeline only ever consumes a real .svg (see dh-book-draft/
  # CLAUDE.md), so this is purely an additional OUTPUT format this
  # function now also accepts -- any caller that still passes a .svg
  # OUTPUT gets the exact same vector file as before.
  if (grepl("\\.png$", output, ignore.case = TRUE)) {
    png(output, width = width_in, height = panel_h_in, units = "in",
        res = 150, type = "cairo", family = getOption("kld.plot.font_family", "Hiragino Mincho ProN"))
  } else {
    svg(output, width = width_in, height = panel_h_in,
        family = getOption("kld.plot.font_family", "Hiragino Mincho ProN"))
  }
  # Bottom/left margins grew from 3.2/3.5 (2026-10-16) to fit the two
  # new rows .draw_axis_direction_arrows() draws between the tick
  # labels and the (now correspondingly pushed-out, see .draw_axes's
  # own title_line) axis title -- a drawn double-arrow row, then an
  # end-text-label row below it. Both ends of this change read the
  # same canvas_scale, so the two stay proportioned together at any
  # canvas_size.
  par(mar = c(5.0, 5.0, 1, 1) * canvas_scale, mgp = c(2, 0.6, 0) * canvas_scale)
  par(cex.axis = canvas_scale, cex.lab = canvas_scale)
  list(canvas_scale = canvas_scale, panel_h_in = panel_h_in)
}

filter_by_edge_strength <- function(d, threshold) {
  # の等虚词类pair跟任何词搭配算出来的association都天然很弱(边际份额太高，
  # "这个词出现"这件事本身没什么信息量)，与keyness/frequency这两个轴衡量
  # 的pair-vs-gender关联无关，因此按pair自身两条边(assoc_w1_to_w2/
  # assoc_w2_to_w1)是否都严格大于阈值过滤，不看轴上的数值；word版表没有
  # 这两列，no-op
  if (is.null(threshold) || is.null(d$assoc_w1_to_w2)) return(d)
  keep <- !is.na(d$assoc_w1_to_w2) & !is.na(d$assoc_w2_to_w1) &
    d$assoc_w1_to_w2 > threshold & d$assoc_w2_to_w1 > threshold
  d[keep, ]
}

filter_by_dispersion <- function(d, threshold) {
  # dispersion(含NA)过低的点/pair只集中在一两首歌里，即使keyness很高也
  # 可能只是偶然；word/pair两版表都有dispersion列
  if (is.null(threshold) || is.null(d$dispersion)) return(d)
  keep <- !is.na(d$dispersion) & d$dispersion > threshold
  d[keep, ]
}

filter_by_node_dispersion <- function(d, threshold) {
  # pair专用，跟filter_by_dispersion互补：那个看pair整体(两词共现)的分布，
  # 这个看两侧entry各自单独的dispersion——pair整体共现分布均匀不代表其中
  # 一侧entry本身不是全语料集中在寥寥几首诗里。两侧都要求>阈值，
  # word版表没有dispersion1/dispersion2，no-op
  if (is.null(threshold) || is.null(d$dispersion1) || is.null(d$dispersion2)) return(d)
  keep <- !is.na(d$dispersion1) & !is.na(d$dispersion2) &
    d$dispersion1 > threshold & d$dispersion2 > threshold
  d[keep, ]
}

filter_by_frequency <- function(d, threshold) {
  # dispersionとは別の判断基準——dispersionは「分布が1-2首に偏っている」
  # を捉える、これはそもそも「出現回数が少なすぎて信頼できない」を捉える
  # (2026-10-15、ユーザーと確認済み：両方残す、片方をもう片方に差し替える
  # のではない)。word/pair両版表ともtotal_freq列を持つのでpair限定の
  # ガードは不要
  if (is.null(threshold) || is.null(d$total_freq)) return(d)
  keep <- !is.na(d$total_freq) & d$total_freq > threshold
  d[keep, ]
}

filter_by_node_frequency <- function(d, threshold) {
  # pair only, mirrors filter_by_node_dispersion: filter_by_frequency
  # looks at the pair's combined total_freq, this looks at each side
  # entry's own freq1/freq2 -- a frequent pair can still have one rare
  # side word. Both sides must exceed threshold. word table has no
  # freq1/freq2, no-op
  if (is.null(threshold) || is.null(d$freq1) || is.null(d$freq2)) return(d)
  keep <- !is.na(d$freq1) & !is.na(d$freq2) &
    d$freq1 > threshold & d$freq2 > threshold
  d[keep, ]
}

.draw_axes <- function(d, ranges, nofreq, axis_label, point_label, x_break = NULL, y_break = NULL,
                        axis_col = "black", kld_variant = NULL) {
  xr <- ranges$xr; yr <- ranges$yr
  plot.new()
  plot.window(xlim = c(xr[1] - ranges$xpad, xr[2] + ranges$xpad),
              ylim = c(yr[1] - ranges$ypad, yr[2] + ranges$ypad))

  kld_sub <- if (nofreq) "nofreq" else "signed"
  # Read-only (never par(mgp = ...) here) -- .draw_axes() runs twice
  # per device in inset mode (source panel then zoomed panel, see
  # .render_with_inset), sharing one .open_canvas() call; mutating
  # par("mgp") here would compound on the second call instead of
  # giving both panels the same offset. title_line pushes the axis
  # title out past where .draw_axis_direction_arrows() now puts its own
  # arrow AND that arrow's own end-text-label row (between the tick
  # labels and the title) -- both read the same base_mgp so they
  # always stay in the same relative order regardless of canvas_scale.
  base_mgp <- par("mgp")
  title_line <- base_mgp[1] + 1.4
  if (nofreq) {
    # 已经是[0,1]量表，不需要再按量级缩放刻度；有断轴的那根轴用
    # .axis_ticks_with_break算出来的不均匀刻度，没有断轴的照常pretty()
    if (!is.null(x_break)) {
      xt <- .axis_ticks_with_break(x_break)
      axis(1, at = xt$at, labels = xt$labels, col = axis_col, col.axis = axis_col, col.ticks = axis_col)
    } else {
      axis(1, at = pretty(xr), col = axis_col, col.axis = axis_col, col.ticks = axis_col)
    }
    if (!is.null(y_break)) {
      yt <- .axis_ticks_with_break(y_break)
      axis(2, at = yt$at, labels = yt$labels, col = axis_col, col.axis = axis_col, col.ticks = axis_col)
    } else {
      axis(2, at = pretty(yr), col = axis_col, col.axis = axis_col, col.ticks = axis_col)
    }
    title(xlab = bquote(plain(KLD)[.(kld_sub)] * "(" * .(point_label) * "→" * .(axis_label) * ")"),
          ylab = bquote(plain(dispersion) * "(" * .(point_label) * ")"), line = title_line)
    usr <- par("usr")
    if (!is.null(x_break)) .draw_axis_break_mark(x_break$main_frac, side = 1, usr = usr)
    if (!is.null(y_break)) .draw_axis_break_mark(y_break$main_frac, side = 2, usr = usr)
  } else {
    # x轴(frequency, TARGET→POINT方向)量级通常远小于y轴(keyness, POINT→TARGET)，
    # 刻度位置不变，显示数字除以x_scale，缩放比例写进标题——用底数^指数
    # (10^.(x_exp))的plotmath写法显示成"×10⁻³"这种科学计数法样式，不是
    # 直接打印0.001这种十进制小数(2026-09-22改，用户反馈原来的纯小数
    # 不够直观)。Only "obs" (unbounded [0,+inf)) actually needs this --
    # "norm" is already rescaled into roughly (-1,1), so the same log10
    # floor() trick just produces a spurious "10^-1" suffix for no
    # readability gain (confirmed with the user, 2026-10-15); KLD_VARIANT
    # being unset (NULL, an older/other caller) keeps the old exponent
    # behavior rather than guessing.
    needs_exponent <- is.null(kld_variant) || identical(kld_variant, "obs")
    x_ticks <- pretty(xr)
    if (needs_exponent) {
      x_exp <- floor(log10(max(abs(xr))))
      x_scale <- 10^x_exp
      axis(1, at = x_ticks, labels = format(x_ticks / x_scale, trim = TRUE),
           col = axis_col, col.axis = axis_col, col.ticks = axis_col)
      axis(2, at = pretty(yr), col = axis_col, col.axis = axis_col, col.ticks = axis_col)
      title(xlab = bquote(plain(KLD)[.(kld_sub)] * "(" * .(point_label) * "→" * .(axis_label) * ")"
                           ~ (10^.(x_exp))),
            ylab = bquote(plain(dispersion) * "(" * .(point_label) * ")"), line = title_line)
    } else {
      axis(1, at = x_ticks, labels = format(x_ticks, trim = TRUE),
           col = axis_col, col.axis = axis_col, col.ticks = axis_col)
      axis(2, at = pretty(yr), col = axis_col, col.axis = axis_col, col.ticks = axis_col)
      title(xlab = bquote(plain(KLD)[.(kld_sub)] * "(" * .(point_label) * "→" * .(axis_label) * ")"),
            ylab = bquote(plain(dispersion) * "(" * .(point_label) * ")"), line = title_line)
    }
  }
  # neutral(理论中性点)如果落在被断轴压缩的那根轴上，要用同一个fwd()换算
  # 成压缩后的位置，不然虚线会画在跟数据脱节的地方
  neutral_x <- if (!is.null(x_break)) x_break$fwd(ranges$neutral) else ranges$neutral
  neutral_y <- if (!is.null(y_break)) y_break$fwd(ranges$neutral) else ranges$neutral
  abline(h = neutral_y, v = neutral_x, col = "#CCCCCC", lty = 2)
  # 中位数是"这批样本对半分"的经验线，跟neutral(population baseline无偏离的
  # 理论中性点)是两个不同概念，不跨图可比——单独画一条实线，不替换neutral的
  # 虚线。d$keyness/d$frequency如果被断轴压缩过，这里读到的已经是压缩后的
  # 坐标(.apply_nofreq_breaks原地替换过)，不需要再额外变换
  dispersion_median <- median(d$keyness)
  dispersion_mean <- mean(d$keyness)
  abline(h = dispersion_median, v = median(d$frequency), col = "#AAAAAA", lty = 1)
  # Mean line + text labels for both (2026-10-15, user-requested) --
  # dispersion's own median/mean only (the X-axis median line above
  # stays unlabeled, out of scope here). Placed at xr[2] (the data
  # range's own right edge, clipped to the panel, xpd=FALSE) rather
  # than the canvas edge -- the color legend/size legend already own
  # the margin past the panel (see .draw_legend/.draw_size_legend),
  # text out there would collide with them.
  abline(h = dispersion_mean, col = "#888888", lty = 3)
  text(xr[2], dispersion_median, labels = sprintf("median %.3f", dispersion_median),
       col = "#AAAAAA", cex = 0.55, adj = c(1, -0.3))
  text(xr[2], dispersion_mean, labels = sprintf("mean %.3f", dispersion_mean),
       col = "#888888", cex = 0.55, adj = c(1, -0.3))
}

# 点の大きさはtotal_freq(n_target+n_not_target)ではなくn_target(target
# 語料Tの中だけでの生頻度)基準——2026-10-17、ユーザー指摘：生頻度を線形
# のままnode sizeに反映すると、低頻度語同士がほぼ同じ大きさになってしまう
# (freq/max_freqの比がほとんどの点で0に近いまま)。対数(底2、+1でfreq=0
# もlog(1)=0に落ちるようにする)を取ってからmax-min正規化(Tに一度も出ない
# 語のlog=0が既にこの正規化の下限なので、実質ただのL(w)/max_v L(v)に
# 簡約できる、ユーザー自身の提示した式通り)。
#
# この正規化計算(log2 + 全体のmax)自体はscripts/patches/visprepare.py
# 側(表生成時、_point_row/_max_log_n_target)で一度だけ行い、結果を
# node_size列としてtsvに直接書き出す——ここ(Rの描画側)では一切計算せず、
# 読み込んだnode_size列をそのままsqrt(面積比例)の式に差し込むだけ
# (2026-10-17、ユーザー確認済み：「在生成表格的部分进行修改避免反复
# 计算」——Shiny appは「描画」を押すたびにこの関数を数千行分呼び直す
# ので、全体のmaxを取る集約計算を毎回Rでやり直すのではなく、表を作る
# 時に1回だけ計算して列に焼き込んでおく)。0.7/0.3/1.1は2026-09-26
# 確認済みの大きさ係数、変更なし。
.size_scale <- function(node_size) 0.7 * (0.3 + 1.1 * sqrt(node_size))

# 扩大面板里点的放大倍数(见render_zoomed_scatter_panel)——取景框越窄，
# 原始倍数(sqrt(zoom_x*zoom_y))越容易冲到20~30倍，封顶在4倍：正常放大
# (框没缩得特别紧)时倍数远低于4，跟原来一样按实际放大程度显示；框缩得
# 特别紧时不再跟着无限变大，最多到4倍(2026-09-26用户反馈确认这个上限，
# 讨论过封顶的是整个面板共用的这一个数，不是每个点各自的最终大小，
# 所以不会出现同一面板里部分点封顶、部分没封顶的不一致)。
.capped_zoom_factor <- function(zoom_factor, cap = 4) min(zoom_factor, cap)

# entry/concept/topic各自的pch(形状)——跟.render_pair_network_svg里network
# node的形状映射(entry=ellipse/concept=box/topic=diamond)对应的散点版：
# entry用空心圆(pch=1)，concept用空心正方形(pch=0，2026-09-22"不是entry的
# field的点改为正方形")，topic用空心菱形(pch=5，呼应network里的diamond)。
# pch<15全部是"只画边框、不填色"的空心符号(2026-09-23"全部改成空心点"，
# 取代之前16/15/18那套实心画法——点数量一多时实心圆互相压住看不清各自的
# dispersion深浅，空心只画边框，重叠时仍能看出轮廓)。只有合并了多种非pair
# field(见draw_point_vs_target()里d$field那段)、每行有自己的field值时才
# 逐点套用；单一point_field(没有d$field列，老缓存tsv/未合并场景)继续用
# 原来"entry圆、非entry方"整图统一一种形状的逻辑，只是形状本身也换成空心。
field_shape_by_type <- function() c(entry = 1, concept = 0, topic = 5)

# 形状图例(.draw_shape_legend)上field名的显示文字——跟run_meta_config_
# plots.R/tune_pair_thresholds_server.R里point_label用的译法(語/概念/部立)
# 保持一致(2026-09-23"轴名中concept改为概念，topic改为部立"、"legend里的
# label也应该改为日语")，图例跟轴标题、跟图内其它文字统一用日语，不留英文
# 原名。
field_label_by_type <- function() c(entry = "語", concept = "概念", topic = "部立")

# ---- 颜色 ----

# |x|归一化渐变：light_col->dark_col，vals可以带正负号(用绝对值定深浅)
.abs_ramp <- function(vals, abs_range, light_col, dark_col) {
  norm <- if (diff(abs_range) > 0) {
    (abs(vals) - abs_range[1]) / diff(abs_range)
  } else {
    rep(0.5, length(vals))
  }
  ramp <- colorRampPalette(c(light_col, dark_col))(101)
  ramp[pmin(100, pmax(1, round(norm * 100))) + 1]
}

# dispersion_difference(D_T - D_R)是带符号的量，.abs_ramp先abs()再映射，
# 正负取值相同大小时颜色完全一样，图例两端数字也因此永远非负(2026-10-17
# 用户发现"legend里没有负值")。中间试过按分位数(ecdf排名)归一化，但这份
# 数据里有36.7%的点dispersion_difference恰好等于0.0(真实相同值，不是
# 浮点误差，shinkokin实测9060/24582)，分位数对这些完全相同的点必然给
# 同一个颜色，这一大块加上其余分散的点，看起来就是"一块黑一块灰"两种
# 色块(2026-10-17用户反馈)，不是真渐变。最终改回2026-10-17用户要的
# "直接用归一化"：对真实的带符号值(signed_range，不取abs)做最朴素的
# 线性min-max映射，风格上跟最早(改成黄蓝之前)的单一灰阶.abs_ramp一样，
# 区别只是这里不取abs()，所以正负在颜色上仍有区别，图例两端标的也是
# 真实的signed_range(可以是负数)，不是之前abs_range那种永远非负的显示。
.linear_ramp <- function(vals, val_range, lo_col, hi_col) {
  norm <- if (diff(val_range) > 0) (vals - val_range[1]) / diff(val_range) else rep(0.5, length(vals))
  ramp <- colorRampPalette(c(lo_col, hi_col))(101)
  ramp[pmin(100, pmax(0, round(norm * 100))) + 1]
}

# 类别配色候选色带，按highlight$type首次出现的顺序依次领取(默认/未命中
# highlight的点固定用蓝色带，不占用这个轮换列表)
category_color_palettes <- function() list(
  c(light = "#F4B6AC", dark = "#7A0C00"),
  c(light = "#C9E9CE", dark = "#164A22"),
  c(light = "#E5D4F0", dark = "#4A1B6B"),
  c(light = "#FCE8B2", dark = "#8A5A00"),
  c(light = "#B8E2E8", dark = "#0D4F5C")
)
default_color_palette <- function() c(light = "#CBD9F0", dark = "#0B2C6B")
# highlight/highlight_pairs都没给(没有多个种类的点，只有一种)时，用这条
# 黑灰渐变代替上面的蓝色——蓝色是"未命中highlight的默认桶"这个语义色，
# 只有存在别的type需要跟它区分时才有意义；根本没有分类时继续用蓝色会
# 让人误以为这也是某种"type"，改用中性的黑灰色(2026-09-22用户反馈)
neutral_color_palette <- function() c(light = "#D9D9D9", dark = "#1A1A1A")
# 未命中highlight/highlight_pairs的"默认"点用这一条从lo(数值低)到
# hi(数值高)的线性渐变(见.linear_ramp)。配色几经调整：先试过蓝/黄两个
# 色相的发散色标，用户否掉了("不用分两个颜色")；2026-10-17最终确认是
# 单一色相的黑色渐变(浅灰->黑)，不是两种颜色。命中highlight某个type的
# 点仍然用该type自己领到的单色相light/dark+abs幅度(.abs_ramp)，保留
# "各type一个颜色"的分类辨识度，没有改动。
gradient_color_palette <- function() c(lo = "#E5E5E5", hi = "#1A1A1A")

# 颜色深浅的量级是dispersion(长表里的type=="dispersion"，衡量该点出现的
# 诗篇分布有多不均匀，跟keyness/frequency是完全独立的另一种量)，不掺
# keyness/frequency。默认(未命中highlight表的点)用固定的蓝色带；highlight
# 表(term, type两列)里每个type各领一条色带、按命中的点各自的dispersion值
# 渐变，覆盖在默认蓝色之上——后出现的type会覆盖先出现的(通常各type的term
# 不重叠，重叠时以highlight表的行序为准)。highlight只管染色，不涉及是否
# 画badge(badge候选是annotate_words这份独立的单列名单，见.select_badge_rows)
# ——两份名单分开维护是2026-09-13定下来的：badge本身的文字已经是它的视觉
# 信号，不需要再用点的颜色重复标一次，混在一起会导致同一批badge一半红
# 一半绿(取决于是否也命中highlight的某个type)这种"颜色没统一"的假象。
# annotate_point_exclude(2026-10-19追加)同样会从highlight的命中结果里
# 摘掉——highlight按读音(point_stem)匹配，跟annotate_words同源的同音
# 词条(比如目标词"え【得】"连带命中"え【枝】"这种)如果不摘掉，会出现
# "红点在inset放大框外"的不一致(zoom_annotated_only的框只认annotate_
# words筛过exclude的点，highlight原本却不认exclude)。
# 返回d$point_col(每个点最终画出来的颜色，供badge/连线复用)和
# type_ramps(具名list，图例要用)。
.compute_point_colors <- function(d, highlight, highlight_pairs = NULL,
                                   annotate_point_exclude = NULL,
                                   highlight_frontier = FALSE,
                                   frontier_n = 15,
                                   category_palettes = category_color_palettes(),
                                   default_palette = default_color_palette(),
                                   neutral_palette = neutral_color_palette(),
                                   gradient_palette = gradient_color_palette()) {
  # Shade now tracks dispersion_difference (D_T - D_R, signed) instead
  # of the pooled dispersion already shown on the Y axis. A tsv
  # regenerated after this column existed but without the aux stats
  # wired for it has the column present but blank ("" -> NA on read);
  # a tsv never regenerated at all (anything but entry-vs-anthology-
  # shinkokin as of 2026-10-15) doesn't have the column AT ALL -- d$
  # dispersion_difference is NULL, not a same-length NA vector, and
  # ifelse(is.na(NULL), 0, NULL) silently returns a ZERO-length vector,
  # not nrow(d) zeroes (confirmed 2026-10-15: broke entry-vs-gender-
  # first4.tsv with "replacement has 0 rows, data has 1655"). Both
  # cases -> uniform shade (0 everywhere), not an error.
  dispersion <- if (is.null(d$dispersion_difference)) {
    rep(0, nrow(d))
  } else {
    ifelse(is.na(d$dispersion_difference), 0, d$dispersion_difference)
  }
  # signed_range is both the real min/max for the legend's two end
  # labels AND now (2026-10-17, reverted off quantile ranking back to
  # plain linear normalization -- see .linear_ramp's own comment) the
  # range the color mapping itself is scaled against. abs_range stays
  # magnitude-only, kept for the per-type highlight ramps below
  # (.abs_ramp), which still encode "how different" via a single hue's
  # light/dark depth, not sign.
  signed_range <- range(dispersion)
  abs_range <- range(abs(dispersion))
  has_categories <- (!is.null(highlight) && nrow(highlight) > 0) ||
    (!is.null(highlight_pairs) && nrow(highlight_pairs) > 0) ||
    isTRUE(highlight_frontier)
  default_palette <- if (has_categories) default_palette else neutral_palette
  # alpha=0.35(其他颜色方案沿用的默认值)在白底上合成灰阶时压缩得特别
  # 厉害：E5E5E5合成后≈246(几乎看不见)，1A1A1A合成后也只有≈175(仍是
  # 浅灰)，整条灰阶的可见反差被压缩到很窄的区间——蓝/黄这种有色相的
  # 配色在同样的alpha下还能凭色相认出来，纯灰阶就不行。这里单独把
  # alpha调高，让灰阶本身的深浅在未叠加时也能看出来(2026-10-17用户反馈
  # "完全不是渐变色，是一块黑一块灰色")。
  point_col <- adjustcolor(.linear_ramp(dispersion, signed_range, gradient_palette["lo"],
                                         gradient_palette["hi"]),
                            alpha.f = 0.85)
  # 图例那条竖条也要用同一套线性映射才能跟正文点的颜色对应得上——这里
  # 直接把60个色块的颜色算好随返回值一起带出去，.draw_legend自己不用再
  # 认识.linear_ramp这些细节，只管拿着一串现成颜色画矩形(跟.draw_legend
  # 之前"自己重新生成一遍渐变色"的写法比，职责更干净，也避免两处对
  # "怎么归一化"各写一份、容易改一边忘另一边)。
  legend_n_steps <- 60
  legend_step_vals <- seq(signed_range[1], signed_range[2], length.out = legend_n_steps)
  legend_cols <- .linear_ramp(legend_step_vals, signed_range, gradient_palette["lo"],
                               gradient_palette["hi"])
  type_ramps <- list()
  if (!is.null(highlight) && nrow(highlight) > 0) {
    point_stem <- sub("【.*$", "", d$point)
    types <- unique(highlight$type)
    for (i in seq_along(types)) {
      type_name <- types[i]
      pal <- category_palettes[[((i - 1) %% length(category_palettes)) + 1]]
      terms <- highlight$term[highlight$type == type_name]
      is_match <- point_stem %in% terms
      if (!is.null(annotate_point_exclude)) {
        is_match <- is_match & !(d$point %in% annotate_point_exclude)
      }
      if (any(is_match)) {
        point_col[is_match] <- .abs_ramp(dispersion[is_match], abs_range, pal["light"], pal["dark"])
      }
      type_ramps[[type_name]] <- pal
    }
  }
  # pair版highlight(hl，比如kondo2011内部拆解出的成分对)：跟highlight同样
  # "只染色，不涉及是否画badge"的定位(2026-09-18确认)，命中判断用d$w1/
  # d$w2这两列(extract_point_vs_target.py对pair_field直接按"|"拆出来的原始
  # 两侧，不是显示用的"・"拼接串)，各自去【...】取stem，再跟highlight_pairs
  # (w1,w2,type三列)做无序集合匹配——不能对d$point按"・"再拆一次：entry
  # 自己的value可能内含"・"(比如"あかし【明し・赤し】"这种多个读音/汉字
  # 并记的写法)，对显示字符串按"・"切会把这种entry自己内部的"・"也切开，
  # 切出三段而不是两段(2026-09-19发现的真实bug，w1/w2这两列就是为了避免
  # 这个问题存在的，不要用point反推)。凡例の名前はCLI引数じゃなくtype列
  # (テーブル側)から取る(2026-09-20確認"hl需要label名，从表格中得到，不是
  # 单纯hl")——highlight(単語版)のtype列と同じ扱い。type列が無い古い名単
  # (2列だけのTSV)にも後方互換で対応：その場合は単一カテゴリ"hl"として
  # 扱う。パレット色は単語版highlightで使った分の続きから割り当てる
  # (両方同時に使っても色が衝突しない)。
  if (!is.null(highlight_pairs) && nrow(highlight_pairs) > 0) {
    point_keys <- mapply(function(a, b) {
      paste(sort(sub("【.*$", "", c(a, b))), collapse = "\t")
    }, d$w1, d$w2)
    pair_types <- if (!is.null(highlight_pairs$type)) unique(highlight_pairs$type) else "hl"
    base_idx <- length(type_ramps)
    for (i in seq_along(pair_types)) {
      type_name <- pair_types[i]
      pal <- category_palettes[[((base_idx + i - 1) %% length(category_palettes)) + 1]]
      rows_i <- if (!is.null(highlight_pairs$type)) {
        highlight_pairs[highlight_pairs$type == type_name, ]
      } else {
        highlight_pairs
      }
      pair_keys <- unique(apply(rows_i[, c("w1", "w2")], 1,
                                 function(r) paste(sort(trimws(r)), collapse = "\t")))
      is_match <- point_keys %in% pair_keys
      if (any(is_match)) {
        point_col[is_match] <- .abs_ramp(dispersion[is_match], abs_range, pal["light"], pal["dark"])
      }
      type_ramps[[type_name]] <- pal
    }
  }
  # highlight_frontier(2026-10-19新增，用户要求"不用highlight_file、用
  # frontier_n和filter配合"之后，又想要"frontier区域也能染色"：不依赖
  # 外部词表，直接取frequency(横轴KLD)从高到低排前frontier_n个，整体染成
  # 一个单独的type("frontier")——跟badge_selection是否真的等于"frontier"
  # 无关，可以单独开，也可以跟words/zoom_box模式的badge并用(颜色只管
  # 染色，不影响选哪些点做badge，跟highlight/highlight_pairs同一个设计
  # 原则)。不设下限/中点(2026-10-19用户明确要求"不要定下限"，
  # frequency_threshold把候选池筛小之后，中点以上的点数会跟着缩水、
  # 凑不够frontier_n个)，直接跟.select_badge_rows的frontier分支同一个
  # "排序后截断前N个"算法，保证只要d里有足够多的点就一定凑够frontier_n。
  if (isTRUE(highlight_frontier)) {
    ord <- order(-d$frequency)
    top_n <- min(frontier_n, nrow(d))
    top_idx <- ord[seq_len(top_n)]
    is_match <- seq_len(nrow(d)) %in% top_idx
    if (!is.null(annotate_point_exclude)) {
      is_match <- is_match & !(d$point %in% annotate_point_exclude)
    }
    base_idx <- length(type_ramps)
    pal <- category_palettes[[(base_idx %% length(category_palettes)) + 1]]
    if (any(is_match)) {
      point_col[is_match] <- .abs_ramp(dispersion[is_match], abs_range, pal["light"], pal["dark"])
    }
    type_ramps[["frontier"]] <- pal
  }
  list(point_col = point_col, signed_range = signed_range, legend_cols = legend_cols,
       type_ramps = type_ramps, default_palette = default_palette)
}

# ---- badge候选筛选 ----

# 四种badge候选筛选方式，靠badge_selection切换：
#   "frontier" -- 按横轴(frequency，轴交换后其实是point->target的signed
#     KLD)从高到低排序，直接取前frontier_n个，不依赖任何外部词表(2026-
#     09-14确认放弃帕累托前沿的方案；2026-09-14到2026-10-19之间曾经是
#     "纵轴dispersion中点以上"，后改"按横轴KLD排序"；2026-10-19又去掉了
#     中点下限这一步，用户明确要求"不要定下限"——配合frequency_threshold
#     这类filter先筛小候选池时，中点以上的点数会跟着候选池一起缩水，
#     经常凑不够frontier_n个；不设下限、直接排序截断，只要候选池本身够
#     大就一定能凑够)。
#   "pairs" -- 仅pair表(有w1/w2列)可用，annotate_pairs/annotate_pairs_file
#     精确指定"就是这几对"(w1,w2两列的名单，无序集合匹配，同annotate_words
#     一样是bare stem，不带【...】)，不是"含有某个词的所有pair"这种模糊
#     匹配——annotate_words对pair表只能框到w1这一侧(见.select_badge_rows
#     "words"分支下面的说明)，想精确点名一个/几个具体pair时用这个模式
#     (2026-09-25用户提出"如果我只想标注一个pair怎么做")。
#   "words"(默认) -- annotate_words有传的话(单列的字符向量，不带类别——
#     badge本身的文字已经是视觉信号，不需要类别配色)，只标这份名单命中的
#     点(不是自动挑出来的极值，而是外部认定的词)；annotate_nontarget_only/
#     annotate_point_exclude是额外的收窄/排除规则(比如"只保留从不出现在
#     target一侧的"、"手动排除同音异字误伤")。对pair表来说，匹配键是
#     sub("【.*$","",d$point)，而pair的point是"w1【…】・w2【…】"，这个
#     regex只切到第一个"【"为止——等于只按w1的读音过滤，同一个w1配的所有
#     w2都会一起命中，不能唯一定位到某一个具体pair(这不是bug，是这个模式
#     本来的设计边界；想精确到单个pair用上面的"pairs"模式)。
#   annotate_words没传时(且不是frontier/pairs模式)退回默认行为：按keyness
#   排序取最极端的label_top_n个。hi/lo分组含义不变：用neutral而不是排名
#   决定某个匹配词该算hi还是lo，只是分组不对应固定颜色(见.badge_color_from_point)。
# pair表里同一对entry两个方向各存一行(w1・w2、w2・w1，坐标/keyness/
# frequency完全相同，见extract_point_vs_target.py)，badge候选按这个键
# 去重，保留先出现的那一行——不然同一对entry会被画成两个内容完全一样的
# badge(用户2026-09-19指出"同两个node的两个方向不必重复绘制")。单词版
# rows没有w1/w2列，原样返回。
.dedup_pair_rows <- function(rows) {
  if (is.null(rows$w1)) return(rows)
  # 极少数pair(实测127,250对里有136对)缺association数据(见
  # pair-entry-entry-association-kokin.tsv构建时的"missing association
  # data"提示——entry×entry这两侧互相cross-compare的cell在array_compute
  # 阶段偶尔会落进errors而不进最终矩阵，具体原因跟pair本身合不合法无关，
  # 属于已知的稀疏例外)，没有v12/v21就画不出方向配色的node+edge小图，
  # 提前滤掉，不要留到graphviz渲染node+edge小图那步才报错(2026-09-19实测：候选一多
  # (zoom_box模式几千个)撞上缺失值的概率显著变大，之前只测15个frontier
  # 候选一直没撞上，不代表这个例外不存在)。
  has_assoc <- !is.na(rows$assoc_w1_to_w2) & !is.na(rows$assoc_w2_to_w1)
  rows <- rows[has_assoc, ]
  key <- paste(pmin(rows$w1, rows$w2), pmax(rows$w1, rows$w2), sep = "\t")
  rows[!duplicated(key), ]
}

# pair行(有w1/w2列)按方向无关的组合键识别，跟.dedup_pair_rows同一个理由：
# 同一对entry两个方向的point字符串不同(w1・w2 vs w2・w1)，不能直接用point
# 列本身当身份；单词/概念/部立行(没有w1/w2列)本来就每行一个point，直接用
# point列。
.point_identity_key <- function(rows) {
  if (!is.null(rows$w1)) paste(pmin(rows$w1, rows$w2), pmax(rows$w1, rows$w2), sep = "\t")
  else rows$point
}

# 散点面板"只画被标注的点"共用的筛选原语——ROWS按是否命中BADGE_ROWS过滤，
# 两处调用点(zoom_annotated_only/repel_annotated_only)各自决定何时调用它，
# 这个函数本身不关心inset/repel这些上层开关。
filter_to_annotated_points <- function(rows, badge_rows) {
  rows[.point_identity_key(rows) %in% .point_identity_key(badge_rows), ]
}

# pair(w1,w2)的无序stem匹配键——去掉各自的【...】取bare stem再排序拼接，
# annotate_pairs/highlight_pairs这类"w1,w2两列的名单"都用bare stem书写
# (跟annotate_words单列名单同一套约定)，不是d$point里"読み【漢字】"那种
# 带括号的显示形式，所以两侧都要sub()一次才能比对上
.pair_stem_key <- function(w1, w2) mapply(function(a, b) {
  paste(sort(sub("【.*$", "", c(a, b))), collapse = "\t")
}, w1, w2)

.select_badge_rows <- function(d, annotate_words, label_top_n, neutral,
                                annotate_nontarget_only, annotate_point_exclude,
                                badge_selection = "words", frontier_n = 15,
                                zoom_box = NULL, annotate_pairs = NULL) {
  if (identical(badge_selection, "frontier")) {
    matched <- .dedup_pair_rows(d)
    matched <- matched[order(-matched$frequency), ]
    if (nrow(matched) > frontier_n) {
      matched <- matched[seq_len(frontier_n), ]
    }
    return(matched)
  }
  if (identical(badge_selection, "zoom_box")) {
    # 扩大范围手动指定时(zoom_box参数)，范围内的点全部annotate，不再按
    # frontier那样只挑keyness最高的frontier_n个——"扩大的范围内点都需要
    # annotate"(用户2026-09-19明确要求)，取景范围本身就是选取条件
    stopifnot(!is.null(zoom_box))
    in_box <- d$frequency > zoom_box$x1 & d$frequency < zoom_box$x2 &
      d$keyness > zoom_box$y1 & d$keyness < zoom_box$y2
    return(.dedup_pair_rows(d[in_box, ]))
  }
  if (identical(badge_selection, "pairs")) {
    # 只对pair表(有w1/w2列)有意义，精确点名"就是这几对"——跟"words"模式
    # 对pair表只能框到w1一侧不同，这里两侧都要匹配上才算命中
    stopifnot(!is.null(d$w1), !is.null(annotate_pairs))
    target_keys <- .pair_stem_key(annotate_pairs$w1, annotate_pairs$w2)
    matched <- d[.pair_stem_key(d$w1, d$w2) %in% target_keys, ]
    return(.dedup_pair_rows(matched))
  }
  if (!is.null(annotate_words)) {
    point_stem <- sub("【.*$", "", d$point)
    matched <- d[point_stem %in% annotate_words, ]
    if (annotate_nontarget_only) matched <- matched[matched$n_target == 0, ]
    # 读音相同、汉字不同的同音异字会被词根匹配一起撞进来，只能靠明确排除
    if (!is.null(annotate_point_exclude)) {
      matched <- matched[!(matched$point %in% annotate_point_exclude), ]
    }
    hi <- matched[matched$keyness > neutral, ]
    lo <- matched[matched$keyness <= neutral, ]
  } else {
    # Ranks by d$frequency, not d$keyness -- since the 2026-10-15 axis
    # swap in draw_point_vs_target() (d$frequency <- d$keyness;
    # d$keyness <- d$dispersion, always applied before this function
    # ever runs), d$frequency is the one still holding the original
    # forward KLD (point->target), while d$keyness now holds
    # dispersion. Default labels pick the forward-KLD extremes (most
    # positive/negative point->target divergence), not dispersion
    # extremes, confirmed with the user 2026-10-15.
    hi <- d[order(-d$frequency), ][seq_len(min(label_top_n, nrow(d))), ]
    lo <- d[order(d$frequency), ][seq_len(min(label_top_n, nrow(d))), ]
  }
  # hi/lo合并成一批一起排版——分开各跑一次的话，两组各自内部不重叠，但
  # hi跟lo之间互相没有约束，仍然可能撞在一起
  .dedup_pair_rows(rbind(hi, lo))
}

# keyness/frequency是(n_target, n_not_target)的确定性函数(同一批freq_mode/
# nofreq/anthology下)——badge候选一多，经常有好几个entry/concept/pair的
# 女/非女计数完全一样(比如全部是"女:1/¬女:0"这种hapax)，落在同一个locator
# 上，会画出好几个内容(计数行)完全雷同、互相挤在一起的badge(用户2026-09-21
# 指出"一个点（locator）对应的condition过多")。这里按(n_target, n_not_target)
# 分组，同一组只留一个代表行(跟.dedup_pair_rows"组内挑一行代表"同一个既有
# 约定，取第一行的point_col/dispersion等列，颜色/keyness/frequency不用另外
# 处理，组内本来就相同)，额外挂一个group_size(组内几个)和members(word模式是
# 字符向量、pair模式是完整子表，供.badge_geometry决定要不要切换成列表/网络
# glyph)。group_size==1的组行为要跟没有这个聚合步骤时完全一致——不能引入任何
# 视觉变化，这是向后兼容的关键验证点。
.aggregate_locator_groups <- function(rows) {
  if (nrow(rows) == 0) return(rows)
  is_pair <- !is.null(rows$w1)
  # 分组键额外带上"种类"(pair: field1+field2组合，比如entry+entry跟
  # entry+concept是两种；word: field本身，entry/concept/topic)——只在合并
  # 了多种pair/多种非pair field的input时才会有实际差异(单一种类时kind_key
  # 对所有行都一样，分组结果跟只按locator分组完全等价，不影响原有行为)。
  # 不这样做的话，不同种类的pair/condition只要恰好共享同一个(n_target,
  # n_not_target)，就会被硬凑成一个名单/网络图混在一起，读者分不清哪个
  # 成员属于哪种(2026-09-23"badge 不同种类的pair和condition保持不聚合")。
  kind_key <- if (is_pair) {
    if (!is.null(rows$field1)) paste(rows$field1, rows$field2, sep = "|") else ""
  } else {
    if (!is.null(rows$field)) rows$field else ""
  }
  key <- paste(rows$n_target, rows$n_not_target, kind_key, sep = "\t")
  idx_by_key <- split(seq_len(nrow(rows)), key)
  # field1/field2(w1/w2各自是entry/concept/topic哪一种)是2026-09-23新增
  # 的列，只在合并了多种pair(entry-entry/entry-concept/entry-topic)的
  # input时才有意义、且才需要靠它选network node形状——老数据/非合并场景
  # 没有这两列(is.null(rows$field1))就不选，跟原来完全一样不受影响。
  pair_member_cols <- c("w1", "w2", "freq1", "freq2", "assoc_w1_to_w2", "assoc_w2_to_w1")
  if (is_pair && !is.null(rows$field1)) pair_member_cols <- c(pair_member_cols, "field1", "field2")
  out_rows <- lapply(idx_by_key, function(idx) {
    rep_row <- rows[idx[1], , drop = FALSE]
    rep_row$group_size <- length(idx)
    rep_row$members <- if (is_pair) {
      list(rows[idx, pair_member_cols])
    } else {
      list(rows$point[idx])
    }
    rep_row
  })
  do.call(rbind, out_rows)
}

# badge填色框(不是空心边框)+深色填充时文字/边框换白保持对比度这套画法
# 照搬scripts/kld/locator_plot.R的draw_corner_locator()。填色直接照抄该点
# 自己实际显示的颜色(rows$point_col)——点是哪个类别的渐变色，badge就跟着
# 哪个类别；点是默认蓝的，badge也是蓝的，不强求所有badge用同一色相。
# 深浅判断(该用白字还是黑字)按亮度算。
.badge_color_from_point <- function(rows) {
  if (nrow(rows) == 0) {
    return(list(fill = character(0), text_col = character(0), border_col = character(0),
                 word_col = character(0)))
  }
  fill <- rows$point_col
  # point_col可能带alpha(默认蓝点是adjustcolor(...,alpha.f=0.35))，跟白色
  # 背景先做alpha合成，再算感知亮度，才能判断"合成后实际看起来"是深是浅
  rgba <- col2rgb(fill, alpha = TRUE)
  a <- rgba["alpha", ] / 255
  r_eff <- rgba["red", ] * a + 255 * (1 - a)
  g_eff <- rgba["green", ] * a + 255 * (1 - a)
  b_eff <- rgba["blue", ] * a + 255 * (1 - a)
  luminance <- (0.299 * r_eff + 0.587 * g_eff + 0.114 * b_eff) / 255
  dark <- luminance < 0.45
  # word_col(框外、画在空白画布上的词文字颜色)固定用不透明版本的同一个
  # 颜色，不跟着text_col按框内深浅换白——词文字不在填色框里，背景永远是
  # 画布底色，用白色会看不见
  word_col <- rgb(rgba["red", ], rgba["green", ], rgba["blue", ], maxColorValue = 255)
  list(fill = fill, word_col = word_col,
       text_col = ifelse(dark, "white", "#111111"),
       border_col = ifelse(dark, "white", "#555555"))
}

# ---- pair网络图(仅.render_with_inset()"扩大图整个范围只有一个locator"
# 这一种场景用，统一走graphviz渲染) ----

# 2026-09-21用户看过几版R原生手绘(圆点、直边、弧线)和graphviz原生渲染的
# 预览后明确选定graphviz的样式("preview_network2.svg是最符合要求的")。
# 2026-09-22改回：pair badge本身(不管是不是一堆pair共享同一locator)日常
# 就是纯文字("w1・w2"，见.badge_geometry)，不再逐个画网络图——网络图只
# 保留给.render_with_inset()"整个扩大图范围内所有候选点汇总起来只对应
# 一个locator"这一种边缘场景，画一整块面板那么大的图(用户原话"如果enlarge
# 的范围里只有一种女n非女n的点，则整合成一个网络图")。用graphviz(neato
# 布局+自己的-Tsvg渲染)产出椭圆node(文字在内部)+双向弧线+node边框配色，
# 整段svg输出当矢量内容，后续通过.splice_network_glyphs()拼接进主图的svg
# 文件。bgcolor=transparent：graphviz默认会给画布填一个不透明白色背景，
# 整段svg拼进主图后会变成一块挡住背后连线/散点的白板(用户反馈"网络图会
# 遮挡其他的线条，希望保持透视")，声明透明背景后graphviz就不画这块背景
# 了，嵌套svg元素本身默认也是透明的，背后的内容能透出来。
# 每个badge的两行(计数框+词，或计数框+名单列表)合并成一个整体框，一起
# 参与排版；内部上下堆叠的相对位置不变——只是这个整体框会被摆到某个最终
# 位置。2026-09-22改回：pair的badge不再画graphviz网络图(node+edge小图)，
# 下方跟word模式一样是纯文字——group_size==1时是一行"w1・w2"(nakaguro
# 分隔)，group_size>1时是"w1・w2"这样的行的列表，跟word模式的单词/名单
# 列表走同一套kind="word"/"word_list"渲染逻辑，不再有专门的"pair_network"
# kind(网络图只在.render_with_inset()"整个扩大图范围只对应一个locator"
# 这个特殊场景下画一次，见该函数，不再是每个badge各画一份)。hl(kondo2011
# 内部拆解pair)是完全独立的另一件事，只影响.compute_point_colors的染色
# (highlight_pairs)，不影响badge选中哪些点(2026-09-18确认："hl是改变点的
# 颜色，badge和frontier+扩大图配合")。has_groups为FALSE(没跑过
# .aggregate_locator_groups，比如aggregate_locators=FALSE)时group_size
# 恒为1，是向后兼容的关键：不聚合时每行都是自己的badge，逐字节维持原样。
# concept点的point字符串是"label（BG-01-163;BG-01-1630）"这种格式
# (scripts/kld-nofreq/src/kld_nofreq/concept_labels.py的_label_group_key，
# 全角括号+分号连接的BG编号列表)——badge里只需要看得懂的label本身，
# BG编号是内部分类代码，2026-09-22用户要求"badge里去掉bg id的显示"；
# 2026-09-23进一步要求"光看label不知道这个概念里有什么词，加3个例词"。
# 例词数据来自SPARQL(entry->concept的成员关系)，R这边没有直接查询能力，
# 靠scripts/analysis/build_concept_display_labels.py离线生成的小查找表
# (原始复合key -> "label（例词1・例词2・例词3）"这样的完整展示文本)。
# 长表自己的condition1_label/condition2_label列也带这份展示文本，但只
# 覆盖"concept单独出现"的行——pair行的pair Element本身从来不生成自己的
# label(见generate_experiment_vector._apply_concept_display_label，只给
# 内层的concept Element本身赋值，不合并成pair的label)，w1/w2两列因此
# 拿不到，所以这里用同一份独立小查找表(按复合key直接查，不依赖长表的
# 哪一列)统一覆盖单点/pair两种情况。查不到(该key在表里没有，或表文件本身
# 不存在——比如还没跑过那个脚本)时退回原来的"只去掉BG编号"行为，不会比
# 改动前更差。
.concept_display_labels_path <- function() {
  Sys.getenv("CONCEPT_DISPLAY_LABELS", "data/processed/concept-display-labels.tsv")
}

.concept_display_label_map <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      path <- .concept_display_labels_path()
      cache <<- if (file.exists(path)) {
        tbl <- read.table(path, sep = "\t", header = FALSE, quote = "",
                           stringsAsFactors = FALSE, encoding = "UTF-8",
                           col.names = c("key", "display"))
        setNames(tbl$display, tbl$key)
      } else {
        character(0)
      }
    }
    cache
  }
})

# TEI本体でXML NameCharの範囲外の文字(例えば"*"=U+002A、"+"=U+002B)は
# xml:id内で"_x<4桁hex>_"形式にエスケープされたまま、measured_value/w1/
# w2までそのまま素通りしてくる(entry2l2entry_exclusion.pyの既存コメント・
# 実データ"の_x002A_【の】"参照) -- 表示の最後の一歩でだけ元の文字に戻す
# (pair keyを"+"で分割するコードはこの関数を通す前の、まだエスケープされた
# ままの文字列を読む -- 通した後の文字列を再度"+"で分割してはいけない、
# "+"自体がエスケープを解いた結果として出てくることがあるため)。
.unescape_xml_id_chars <- function(s) {
  vapply(s, function(x) {
    if (is.na(x) || !grepl("_x[0-9A-Fa-f]{4}_", x)) return(x)
    matches <- unique(regmatches(x, gregexpr("_x[0-9A-Fa-f]{4}_", x))[[1]])
    for (mt in matches) {
      ch <- intToUtf8(strtoi(substr(mt, 3, 6), base = 16L))
      x <- gsub(mt, ch, x, fixed = TRUE)
    }
    x
  }, character(1), USE.NAMES = FALSE)
}

# 部立大分類(budatemajor、field_label_by_type()の"topic")の表示は【】では
# なく『』で囲む -- 元データのbudatemajorラベル自体には括弧が一切ない
# (例"その他")、entry側の既存の【漢字】記法と視覚的に区別するための表示
#専用の追加(2026-10-18ユーザー要望)。FIELD_TYPEがこの文字列以外(NA含む)
# の時は何もしない -- 複合済みの"w1・w2"文字列をもう一度この関数に通す
# 呼び出し(後述の"word"フォールバック)はFIELD_TYPEにNAを渡し、二重に
# 括弧で囲まれることを防ぐ。
.wrap_topic_label <- function(text, field_type) {
  is_topic <- !is.na(field_type) & field_type == "部立大分類"
  ifelse(is_topic & !is.na(text), paste0("『", text, "』"), text)
}

# pair展示用的w1/w2连接符——budate(部立大分類)pair用"，"，其余(entry+entry/
# entry+concept等)用"+"(2026-10-18用户要求，取代原来统一用"・"的写法)。
# FIELD1/FIELD2两边只要有一个是部立大分類就算budate pair(budate pair里
# 非词汇的那一侧固定是field1，但这里不假设具体是哪一侧，两边都检查)。
.pair_display_sep <- function(field1, field2) {
  is_budate <- (!is.na(field1) & field1 == "部立大分類") | (!is.na(field2) & field2 == "部立大分類")
  ifelse(is_budate, "，", "+")
}

.concept_display_text <- function(s, field_type = NA_character_) {
  hit <- .concept_display_label_map()[s]
  text <- ifelse(!is.na(hit), hit, sub("（BG-[^）]*）", "", s))
  text <- .unescape_xml_id_chars(text)
  text <- .wrap_topic_label(text, field_type)
  text
}

# 标注点表格(app.R「表」tab)展示用——跟图上badge文字同一套清理规则
# (.concept_display_text自身：还原_xHHHH_转义、部立大分類类用『』包起来)，
# 分别套到point(非pair)/w1+w2(pair)这两种给人看的列上。pair行的point列
# 本身是技术性的"+"连接key(.row_id之类下游代码可能还要靠它拆分pair，
# annotated_table是这个函数唯一的消费者，目前没有下游再解析，但保持
# point列原样、只处理w1/w2，跟图上的约定一致——图上badge文字从来也不是
# 拿point列直接显示，是从w1/w2现拼的)。
.format_display_columns <- function(d) {
  if (!is.null(d$w1)) {
    field1 <- if (!is.null(d$field1)) d$field1 else NA_character_
    field2 <- if (!is.null(d$field2)) d$field2 else NA_character_
    d$w1 <- .concept_display_text(d$w1, field1)
    d$w2 <- .concept_display_text(d$w2, field2)
  } else if (!is.null(d$point)) {
    d$point <- .concept_display_text(d$point, d$field)
  }
  d
}

# 合并pair-entry-entry/pair-entry-concept/pair-entry-topic这几种独立的
# enriched tsv到一张图时用(见draw_point_vs_target()的input参数)——不能直接
# rbind：不同批次生成的文件列可能不完全一致(比如field1/field2是2026-09-23
# 才加的列，已经落盘的旧缓存tsv还没有)。取所有表的列名并集，缺的列补NA
# 再rbind，缺field1/field2时网络图节点形状按.render_pair_network_svg里的
# 缺省ellipse处理，不会报错。
.rbind_fill <- function(dfs) {
  all_cols <- unique(unlist(lapply(dfs, names)))
  dfs2 <- lapply(dfs, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    for (col in missing_cols) df[[col]] <- NA
    df[all_cols]
  })
  do.call(rbind, dfs2)
}

.badge_geometry <- function(rows, value_label, not_value_label, badge_cex, fam, label_gap_in,
                              max_group_members = 50) {
  is_pair <- !is.null(rows$w1)
  has_groups <- !is.null(rows$group_size)
  lapply(seq_len(nrow(rows)), function(i) {
    count_label <- sprintf("%s:%d / %s:%d", value_label, rows$n_target[i],
                            not_value_label, rows$n_not_target[i])
    count_w_in <- strwidth(count_label, cex = badge_cex, family = fam, units = "inches")
    count_h_in <- strheight(count_label, cex = badge_cex, family = fam, units = "inches")
    box_w_in <- count_w_in + count_w_in * 0.35
    box_h_in <- count_h_in + count_h_in * 1.2
    group_size <- if (has_groups) rows$group_size[i] else 1

    # pair模式下方文字统一是"w1・w2"(nakaguro分隔)，跟word模式共用word/
    # word_list这两种kind——group_size==1时members_i是只含这一个pair的
    # 单行子表(没跑过.aggregate_locator_groups()时，比如aggregate_locators=
    # FALSE，从rows自己的列现凑一份等价单行子表)，paste出来的members长度
    # 恰好是1，落进下面word分支；group_size>1时members_i是组内所有pair的
    # 子表，paste出来是多行"w1・w2"字符串，走word_list分支——跟word模式的
    # 名单列表(2026-09-21新增的功能)是同一套渲染/截断逻辑，不用为pair单独
    # 再写一套。
    if (is_pair) {
      pair_cols <- intersect(c("w1", "w2", "field1", "field2"), names(rows))
      members_i <- if (has_groups) {
        rows$members[[i]]
      } else {
        rows[i, pair_cols, drop = FALSE]
      }
      field1 <- if (!is.null(members_i$field1)) members_i$field1 else NA_character_
      field2 <- if (!is.null(members_i$field2)) members_i$field2 else NA_character_
      members <- paste0(.concept_display_text(members_i$w1, field1),
                         .pair_display_sep(field1, field2),
                         .concept_display_text(members_i$w2, field2))
    } else if (group_size > 1) {
      members <- .concept_display_text(rows$members[[i]], rows$field[i])
    } else {
      members <- NULL
    }

    if (!is.null(members) && length(members) > 1) {
      shown <- if (length(members) > max_group_members) {
        c(members[seq_len(max_group_members)],
          sprintf("…計%d件", length(members)))
      } else {
        members
      }
      line_w_in <- vapply(shown, strwidth, numeric(1), cex = badge_cex, family = fam, units = "inches")
      line_h_in <- vapply(shown, strheight, numeric(1), cex = badge_cex, family = fam, units = "inches")
      return(list(kind = "word_list", word_list = shown, count_label = count_label,
                   box_w_in = box_w_in, box_h_in = box_h_in,
                   combined_w_in = max(max(line_w_in), box_w_in),
                   combined_h_in = box_h_in + label_gap_in + sum(line_h_in) + label_gap_in * length(shown)))
    }
    # ここに来る時点でis_pairなら(group_size==1でも)membersは常に非NULLで
    # 既に.concept_display_text済み("w1・w2"結合後の文字列) -- 二重に
    # 『』を付けないようFIELD_TYPEはNAを渡す(.wrap_topic_labelはNAなら
    # 素通り)。!is_pairならここに来る時点でgroup_size==1が保証されている
    # (group_size>1は上のword_list分岐がreturn済み)ので、membersはNULL、
    # rows$point[i]を生のまま.field[i]付きで処理する。
    word_field <- if (is_pair) NA_character_ else rows$field[i]
    word <- .concept_display_text(if (!is.null(members)) members[1] else rows$point[i], word_field)
    word_w_in <- strwidth(word, cex = badge_cex, family = fam, units = "inches")
    word_h_in <- strheight(word, cex = badge_cex, family = fam, units = "inches")
    list(kind = "word", word = word, count_label = count_label,
         box_w_in = box_w_in, box_h_in = box_h_in,
         combined_w_in = max(word_w_in, box_w_in),
         combined_h_in = box_h_in + label_gap_in + word_h_in)
  })
}

# ---- badge排版：两种可切换算法 ----

# 力导向repel(ggrepel:::repel_boxes2)——2026-09-12一度被"从画布左上角开始
# 的确定性网格"取代、依赖也从flake.nix移除，2026-09-13又作为badge_layout=
# "repel"选项重新引入，参数沿用移除前最后一版：种子位置摆在点的左上方，
# force_push=1e-4(整体散开)/force_pull=1e-6(适度拉回避免离点太远)。
.layout_badges_repel <- function(combined_w, combined_h, px_in, py_in, pin, badge_gap_in) {
  n <- length(combined_w)
  x1 <- px_in - 0.01 * pin[1] - combined_w
  y1 <- py_in + badge_gap_in
  box_pad_y_in <- badge_gap_in * 1.5
  boxes <- cbind(x1 = x1, y1 = y1 - box_pad_y_in,
                 x2 = x1 + combined_w, y2 = y1 + combined_h + box_pad_y_in)
  repelled <- ggrepel:::repel_boxes2(
    data_points = cbind(px_in, py_in), point_size = rep(0.02, n),
    point_padding_x = badge_gap_in, point_padding_y = badge_gap_in,
    boxes = boxes, xlim = c(0, pin[1]), ylim = c(0, pin[2]),
    hjust = rep(0.5, n), vjust = rep(0.5, n),
    force_push = 1e-4, force_pull = 1e-6, max_time = 3,
    max_overlaps = 20, max_iter = 3000, direction = "both"
  )
  list(x1 = repelled$x, y1 = repelled$y + box_pad_y_in)
}

# 01矩阵里找面积最大的全空(FALSE)矩形——经典"柱状图内最大矩形"逐行算法，
# O(nrow*ncol)，不需要额外依赖。occupied第1行对应画布顶部。
.largest_empty_rect <- function(occupied) {
  nr <- nrow(occupied); nc <- ncol(occupied)
  heights <- integer(nc)
  best <- list(area = 0, r1 = 1, r2 = 0, c1 = 1, c2 = 0)
  for (r in seq_len(nr)) {
    heights <- ifelse(occupied[r, ], 0L, heights + 1L)
    hh <- c(heights, 0L)
    stack <- integer(0)
    for (c in seq_len(nc + 1)) {
      while (length(stack) > 0 && hh[stack[length(stack)]] >= hh[c]) {
        top <- stack[length(stack)]
        stack <- stack[-length(stack)]
        width <- if (length(stack) == 0) c - 1 else c - stack[length(stack)] - 1
        area <- hh[top] * width
        if (area > best$area) {
          best <- list(area = area, r1 = r - hh[top] + 1, r2 = r,
                       c1 = if (length(stack) == 0) 1 else stack[length(stack)] + 1,
                       c2 = c - 1)
        }
      }
      stack <- c(stack, c)
    }
  }
  best
}

# 新默认布局：每个badge各自单独找一块刚好够放下它(w+margin, h+margin)的
# 空白——不预先把画布切成固定几块，而是按row_order(keyness极端程度)顺序，
# 每放一个badge就重新算一次"此刻最大的全空矩形"，在其中挑离这个badge自己
# 的点最近的位置摆放，再把这块地方(含margin)标记为占用，才轮到下一个
# badge重新找(2026-09-14确认"不是4个区块，而是更自由的空白小空间，每个
# 小空间可以放一个badge")。天然贴合每个badge自己的实际尺寸，不会出现
# "预先切出来的区块比某个badge还窄/矮"的问题，也不需要额外的容量估算/
# 区块分配逻辑——固定4区块那版实测下来，无论怎么调整容量估算都还是会有
# badge被分配到偏小的区块里溢出，根源是提前把空白切死了。
.layout_badges_empty_region <- function(combined_w, combined_h, row_order,
                                          px_in, py_in,
                                          all_px_in, all_py_in, point_radius_in,
                                          pin, badge_gap_in, res = 60) {
  n <- length(combined_w)
  cell_w <- pin[1] / res
  cell_h <- pin[2] / res
  occupied <- matrix(FALSE, nrow = res, ncol = res)
  cx_idx <- pmin(res, pmax(1, ceiling(all_px_in / cell_w)))
  cy_idx <- pmin(res, pmax(1, ceiling((pin[2] - all_py_in) / cell_h)))
  pad_x <- max(1, ceiling(point_radius_in / cell_w))
  pad_y <- max(1, ceiling(point_radius_in / cell_h))
  for (k in seq_along(cx_idx)) {
    xr <- max(1, cx_idx[k] - pad_x):min(res, cx_idx[k] + pad_x)
    yr <- max(1, cy_idx[k] - pad_y):min(res, cy_idx[k] + pad_y)
    occupied[yr, xr] <- TRUE
  }
  # 画布右上角是.draw_size_legend/.draw_legend固定占用的区域(图例在badge
  # 之后才画，找空白矩形时不知道那块地方"看起来空"其实已经被图例占了)，
  # 保守地把这块也标记成占用，不然badge可能会被排进去、跟图例文字叠在一起
  # (2026-09-14用4区块方案实测后发现的重叠)
  legend_w_in <- 0.75; legend_h_in <- 2.2
  lx1 <- max(1, ceiling((pin[1] - legend_w_in) / cell_w))
  ly2 <- max(1, ceiling(legend_h_in / cell_h))
  occupied[1:min(res, ly2), lx1:res] <- TRUE

  x1v <- numeric(n); y1v <- numeric(n)
  for (i in row_order) {
    best <- .largest_empty_rect(occupied)
    region <- c((best$c1 - 1) * cell_w, pin[2] - best$r2 * cell_h,
                best$c2 * cell_w, pin[2] - (best$r1 - 1) * cell_h)
    w <- combined_w[i] + badge_gap_in
    h <- combined_h[i] + badge_gap_in
    # 在区域范围内尽量贴近这个badge自己的点摆放(点离得越近连线越短)，
    # 超出区域(badge比这块空白还大，极端情况)时退化成clamp到边界，不做
    # 额外处理
    x1 <- max(region[1], min(px_in[i] - w / 2, region[3] - w))
    y1 <- max(region[2], min(py_in[i] - h / 2, region[4] - h))
    x1v[i] <- x1
    y1v[i] <- y1

    occ_c1 <- max(1, floor(x1 / cell_w) + 1)
    occ_c2 <- min(res, ceiling((x1 + w) / cell_w))
    occ_r1 <- max(1, floor((pin[2] - (y1 + h)) / cell_h) + 1)
    occ_r2 <- min(res, ceiling((pin[2] - y1) / cell_h))
    occupied[occ_r1:occ_r2, occ_c1:occ_c2] <- TRUE
  }
  list(x1 = x1v, y1 = y1v)
}

# ---- badge/连线绘制 ----

# 靠badge的那一段永远轴对齐(纯横或纯竖)，垂直/水平戳进badge外框最近的边；
# 不允许从下边(南侧)接入——点在badge下方本该走南侧入口的情况改走东/西侧。
# 靠原点的那一段允许斜着走，跟第二段之间天然形成钝角(不是90度直角)。
# 全程在inches(等比例)空间算角度，避免x/y轴数量级不同把角度算歪。入口
# 位置固定：东/西侧落在上方1/4处(大致对齐计数框)，南/北侧落在水平正中央。
#
# box_w_in是"有色计数框"自己的宽度，跟combined_w(计数框/下方词语文字两者
# 取较宽的那个，见.badge_geometry)不是一回事——词语文字比计数框宽的时候
# (比如"きりぎりす【蛬・蟋蟀】"这种长词)，combined_w会比计数框本身宽，
# 而计数框是在这个更宽的combined区域里居中摆放的(.draw_badge_box里
# cx <- combined_x1_in + combined_w/2，框以cx为中心画)。东/西侧入口如果
# 直接锚定在combined_x1_in/combined_x2_in(整体区域的左右边缘)，会落在
# 比计数框本身更靠外的地方，看起来像线没接到"有色部分"上，只接到了
# 看不见的、更宽的整体边框(2026-09-15确认的问题："我说的badge是指上方
# 有色部分的边框...但是上方有色部分没有连接")。改成先把计数框自己居中
# 换算出来的左右边缘算出来，东/西侧入口锚定在这个边缘上，跟
# .draw_badge_box实际画色块的位置完全对齐。
.draw_connector_line <- function(px_in, py_in, combined_x1_in, combined_y1_in,
                                  combined_w, combined_h, box_w_in, box_h_in, col,
                                  in_to_data_x, in_to_data_y) {
  combined_x2_in <- combined_x1_in + combined_w
  combined_y2_in <- combined_y1_in + combined_h
  box_cx_in <- (combined_x1_in + combined_x2_in) / 2
  box_cy_in <- (combined_y1_in + combined_y2_in) / 2
  box_x1_in <- box_cx_in - box_w_in / 2
  box_x2_in <- box_cx_in + box_w_in / 2
  enters_vertically <- abs(py_in - box_cy_in) >= abs(px_in - box_cx_in)
  if (enters_vertically && py_in < box_cy_in) enters_vertically <- FALSE
  if (enters_vertically) {
    anchor_x_in <- box_cx_in
    anchor_y_in <- combined_y2_in
  } else {
    # 东/西侧入口要落在计数框本身的上1/4处，不是整个combined box(计数框+
    # 下面的词/node-edge小图)的上1/4——pair版下半部分(node+edge+两个对角
    # label)比单词版的一行文字高很多，combined_h暴涨之后，"整体高度的25%"
    # 会大幅越过计数框、扎进下面的小图区域(用户2026-09-19指出"整体
    # annotation的高度增加导致连线不再在计数框的位置")。改成按box_h_in
    # (只算计数框自己的高度)算，不管下半部分多高，入口位置始终跟计数框
    # 对齐，跟原来的设计意图一致。
    anchor_y_in <- combined_y2_in - 0.25 * box_h_in
    anchor_x_in <- if (px_in <= box_cx_in) box_x1_in else box_x2_in
  }
  dx_in <- anchor_x_in - px_in
  dy_in <- anchor_y_in - py_in
  bend_frac <- 0.6
  if (enters_vertically) {
    mid_x_in <- anchor_x_in
    mid_y_in <- py_in + dy_in * bend_frac
  } else {
    mid_y_in <- anchor_y_in
    mid_x_in <- px_in + dx_in * bend_frac
  }
  lines(in_to_data_x(c(px_in, mid_x_in, anchor_x_in)),
        in_to_data_y(c(py_in, mid_y_in, anchor_y_in)),
        col = col, lwd = 0.8, xpd = NA)
}

.draw_badge_box <- function(g, combined_x1_in, combined_y1_in, combined_w, combined_h,
                             fill, border, text_col, word_col, badge_cex, fam,
                             label_gap_in, x_per_in, y_per_in, in_to_data_x, in_to_data_y) {
  cx_in <- combined_x1_in + combined_w / 2
  cx <- in_to_data_x(cx_in)
  count_top_in <- combined_y1_in + combined_h
  count_bot_in <- count_top_in - g$box_h_in
  y_top <- in_to_data_y(count_top_in)
  box_w <- g$box_w_in * x_per_in
  box_h <- g$box_h_in * y_per_in

  rect(cx - box_w / 2, y_top - box_h, cx + box_w / 2, y_top,
       col = fill, border = border, xpd = NA)
  text(cx, y_top - box_h / 2, labels = g$count_label, cex = badge_cex,
       family = fam, col = text_col, xpd = NA)
  # Word text itself stays fixed black, not word_col (the point's own
  # opaque color -- still used for the connector line below, unchanged)
  # -- badge selection no longer tracks the color ramp's dark end (it's
  # picked by forward-KLD extremes now, see .select_badge_rows), so a
  # labeled word can easily land near abs(dispersion_difference)==0,
  # the ramp's lightest/near-white end, which is unreadable as text.
  # Confirmed with the user 2026-10-15.
  word_text_col <- "#111111"
  if (identical(g$kind, "word_list")) {
    y_in <- count_bot_in - label_gap_in
    for (line in g$word_list) {
      h_in <- strheight(line, cex = badge_cex, family = fam, units = "inches")
      text(cx, in_to_data_y(y_in), labels = line,
           cex = badge_cex, family = fam, col = word_text_col, adj = c(0.5, 1), xpd = NA)
      y_in <- y_in - h_in - label_gap_in
    }
  } else {
    text(cx, in_to_data_y(count_bot_in - label_gap_in), labels = g$word,
         cex = badge_cex, family = fam, col = word_text_col, adj = c(0.5, 1), xpd = NA)
  }
}

# 编排：算geometry -> 挑布局算法 -> 逐个画连线+badge。两种布局算法产出
# 同样格式的位置(list(x1=,y1=))，后续连线/绘制代码不用为两种模式各写一份。
.draw_badges <- function(rows, fill, border, text_col, word_col, badge_layout,
                          d_all, value_label, not_value_label, badge_cex, fam,
                          label_gap_in, x_per_in, y_per_in, data_to_in_x, data_to_in_y,
                          in_to_data_x, in_to_data_y, pin, neutral, max_group_members = 50,
                          tie_boxes = list()) {
  n <- nrow(rows)
  if (n == 0) return(invisible(NULL))
  if (length(fill) == 1) fill <- rep(fill, n)
  if (length(border) == 1) border <- rep(border, n)
  if (length(text_col) == 1) text_col <- rep(text_col, n)
  if (length(word_col) == 1) word_col <- rep(word_col, n)

  geoms <- .badge_geometry(rows, value_label, not_value_label, badge_cex, fam, label_gap_in,
                            max_group_members = max_group_members)
  combined_w <- vapply(geoms, `[[`, numeric(1), "combined_w_in")
  combined_h <- vapply(geoms, `[[`, numeric(1), "combined_h_in")
  box_w <- vapply(geoms, `[[`, numeric(1), "box_w_in")
  box_h <- vapply(geoms, `[[`, numeric(1), "box_h_in")
  px_in <- data_to_in_x(rows$frequency)
  py_in <- data_to_in_y(rows$keyness)
  badge_gap_in <- label_gap_in
  row_order <- order(-abs(rows$keyness - neutral))

  resolved <- if (identical(badge_layout, "repel")) {
    .layout_badges_repel(combined_w, combined_h, px_in, py_in, pin, badge_gap_in)
  } else {
    all_px_in <- data_to_in_x(d_all$frequency)
    all_py_in <- data_to_in_y(d_all$keyness)
    point_radius_in <- 0.04
    .layout_badges_empty_region(combined_w, combined_h, row_order, px_in, py_in,
                                 all_px_in, all_py_in, point_radius_in, pin, badge_gap_in)
  }

  # 先画完全部15个badge的框，连线放在所有框都画完之后统一画——不是每个
  # badge各自"先画线、再画自己的框"这样交替(之前的写法)。原因：badge挤得
  # 很密时，后画的某个badge的框经常会正好压在前一个badge连线的必经路径
  # 上(尤其这批候选点很多是"女:1/¬女:0"这种低频词，KLD算出来的坐标
  # 完全相同，十几条线全从同一个点发散，密集穿过整片badge网格)，把那条
  # 线中间一段盖住，看起来像是"线没接到框上"——线本身其实精确连到了自己
  # 框的边缘，只是被另一个无关的、后画的框挡住了中间一截(2026-09-15确认
  # 的问题，用单独隔离的badge验证过锚点坐标是对的，问题出在图层顺序)。
  # 线全部放到最后画，保证不会被任何框盖住，从起点到终点整条线都看得见；
  # 代价是线的末端会轻微盖在自己框的边框上，但这正是想要的"看得出连上了"
  # 的效果。
  for (i in seq_len(n)) {
    .draw_badge_box(geoms[[i]], resolved$x1[i], resolved$y1[i], combined_w[i], combined_h[i],
                     fill[i], border[i], text_col[i], word_col[i], badge_cex, fam,
                     label_gap_in, x_per_in, y_per_in, in_to_data_x, in_to_data_y)
  }
  # ROWSの各点がtie_boxes(.compute_x_tie_boxesの出力、{key: {x1,y1,x2,
  # y2}})のどれかに属していれば、連結線の出発点はその点自身ではなく箱の
  # 辺にする(2026-10-17、ユーザー確認済み：「線がもし框なら、框のedgeから
  # 出発。もし点なら、そのまま」)。出発点を「箱の中でbadge自身の位置に
  # 一番近い点」にクランプするだけで済む -- badgeは普通箱の外にあるので、
  # x/yそれぞれ独立に[box_x1,box_x2]/[box_y1,box_y2]へクランプした結果は
  # そのまま箱の辺(角を含む)に乗る、標準的な「矩形上の最近点」の求め方。
  start_x_in <- px_in
  start_y_in <- py_in
  if (length(tie_boxes) > 0) {
    keys <- paste(rows$n_target, rows$n_not_target, sep = "/")
    for (i in seq_len(n)) {
      box <- tie_boxes[[keys[i]]]
      if (is.null(box)) next
      box_x_in <- data_to_in_x(c(box$x1, box$x2))
      box_y_in <- data_to_in_y(c(box$y1, box$y2))
      badge_cx_in <- resolved$x1[i] + combined_w[i] / 2
      badge_cy_in <- resolved$y1[i] + combined_h[i] / 2
      start_x_in[i] <- min(max(badge_cx_in, min(box_x_in)), max(box_x_in))
      start_y_in[i] <- min(max(badge_cy_in, min(box_y_in)), max(box_y_in))
    }
  }
  for (i in seq_len(n)) {
    # 连线颜色用word_col(不透明版本)而不是rows$point_col(默认蓝点带
    # alpha=0.35半透明)——连线也是画在空白画布上(道理跟word_col本身的
    # 那条注释一样)，默认类别的半透明蓝配上lwd=0.8的细线，在扩大图这种
    # 点离badge经常隔得比原图更远(empty_region算法找到的空位不一定靠近
    # 点)的场景下，实际渲染出来几乎看不见，容易被误以为"连线没接上"
    # (2026-09-15确认的问题："连线没有接到badge"，根源是可见度不是没画)
    .draw_connector_line(start_x_in[i], start_y_in[i], resolved$x1[i], resolved$y1[i],
                          combined_w[i], combined_h[i], box_w[i], box_h[i], word_col[i],
                          in_to_data_x, in_to_data_y)
  }
  invisible(NULL)
}

# ---- 网络badge的矢量拼接(svg文本级后处理) ----

# 画两个隐形校准标记(用极罕见颜色，人眼看不出来是"多画的东西")，画完之后
# 从最终写盘的svg文件里反查这两个标记实际落在什么原始坐标——因为
# grconvertX/Y(..., to="device"/"ndc")报告的坐标跟cairo实际写进svg文件的
# 坐标对不上(2026-09-21实测：偏差量还跟par(mar=...)设置有关，不是固定
# offset，没找到能直接从R内部拿到这个换算的可靠办法)，只能反过来"画一个
# 已知inches坐标的标记，事后从文件里读它实际被画在哪"这样反解仿射变换。
# 理论上scale应该精确是72(x)/-72(y)(1英寸=72pt，svg的y轴向下、R的inches
# 坐标y轴向上，符号相反)，.splice_network_glyphs()会用两个标记点反解出来
# 的scale校验这个假设，不是就直接报错，不能悄悄套用错的换算继续跑。
.draw_svg_calibration_marks <- function(in_to_data_x, in_to_data_y, pin) {
  mark_col <- "#010203"
  # rect()的坐标要用当前panel的data坐标系，不是inches——标记本身物理上要
  # 很小(0.001英寸见方)，但之前直接把0.001当成data坐标的偏移量硬加/减，
  # 隐含假设了"data坐标数量级接近1"，对多数测试数据(keyness大约0~10)凑巧
  # 成立，但真实数据的x轴(frequency)量级常常只有0.001~0.01(缩放放大图的
  # 可视范围可能窄到只有0.002左右)，同样写死0.001在这种数据上会占掉几乎
  # 整个可视范围，标记从"一个小点"变成"一条贯穿大半个面板的粗线"，进而
  # 让.splice_network_glyphs()反解出来的scale严重偏离72(2026-09-21用真实
  # 数据+改成上下叠放布局的缩放面板才踩出来，不是这次布局改动引入的，是
  # 一直存在、之前凑巧没暴露的bug)。改成用in_to_data_*闭包自己的斜率(1
  # 英寸对应多少个data单位)换算物理偏移量，不管当前data坐标量级多大都成立。
  dx_per_in <- in_to_data_x(1) - in_to_data_x(0)
  dy_per_in <- in_to_data_y(1) - in_to_data_y(0)
  mark_half_in <- 0.001
  marks_in <- list(c(0.15, 0.15), c(max(pin[1] - 0.15, 0.15), max(pin[2] - 0.15, 0.15)))
  for (p in marks_in) {
    x0 <- in_to_data_x(p[1]); y0 <- in_to_data_y(p[2])
    dx <- mark_half_in * dx_per_in; dy <- mark_half_in * dy_per_in
    rect(x0 - dx, y0 - dy, x0 + dx, y0 + dy, col = mark_col, border = NA, xpd = NA)
  }
  marks_in
}

# pending(.draw_badges()收集的待拼接清单，每条含graphviz的svg_txt+目标
# inches区域)拼接进已经写盘的主svg文件——graphviz的svg内容整段当成一个
# 嵌套<svg x y width height viewBox>元素插入主文件(SVG标准支持svg嵌套，
# 嵌套元素会按自己的viewBox自动缩放去填满给定的width/height，不需要手算
# transform矩阵)，取自己的<g ...>到</svg>之间的内容、丢掉graphviz自己的
# <?xml?>/<!DOCTYPE>/外层<svg>标签。marks_in是.draw_svg_calibration_marks()
# 画标记时用的同一组inches坐标，用来反解仿射变换。
.splice_network_glyphs <- function(output, pending, marks_in) {
  if (length(pending) == 0) return(invisible(NULL))
  if (!grepl("\\.svg$", output, ignore.case = TRUE)) {
    stop("badge_layout=\"network\" needs svg output (graphviz glyphs are ",
         "spliced in as SVG XML text) -- png output doesn't support this badge layout")
  }
  svg_txt <- paste(readLines(output, encoding = "UTF-8", warn = FALSE), collapse = "\n")
  lines_all <- strsplit(svg_txt, "\n")[[1]]
  mark_lines <- grep("0.392157%, 0.784314%, 1.176471%", lines_all, fixed = TRUE, value = TRUE)
  if (length(mark_lines) != 2) {
    stop("找不到svg校准标记(应该有2个)，网络badge拼接失败——可能是svg输出格式变了")
  }
  get_center <- function(line) {
    d_str <- sub('.*d="([^"]*)".*', "\\1", line)
    nums <- as.numeric(regmatches(d_str, gregexpr("[0-9]+\\.?[0-9]*", d_str))[[1]])
    c(mean(nums[c(1, 3)]), mean(nums[c(2, 6)]))
  }
  raw1 <- get_center(mark_lines[1]); raw2 <- get_center(mark_lines[2])
  in1 <- marks_in[[1]]; in2 <- marks_in[[2]]
  scale_x <- (raw2[1] - raw1[1]) / (in2[1] - in1[1])
  scale_y <- (raw2[2] - raw1[2]) / (in2[2] - in1[2])
  if (abs(scale_x - 72) > 0.5 || abs(scale_y + 72) > 0.5) {
    stop(sprintf("svg校准标记反解出的scale(%.3f, %.3f)偏离预期的(72, -72)太多，网络badge拼接的坐标换算前提不成立，需要重新核对",
                 scale_x, scale_y))
  }
  offset_x <- raw1[1] - scale_x * in1[1]
  offset_y <- raw1[2] - scale_y * in1[2]
  affine_x <- function(ix) scale_x * ix + offset_x
  affine_y <- function(iy) scale_y * iy + offset_y

  lines_all <- lines_all[!grepl("0.392157%, 0.784314%, 1.176471%", lines_all, fixed = TRUE)]

  blocks <- vapply(pending, function(p) {
    x1_raw <- affine_x(p$x1_in)
    y_top_raw <- affine_y(p$y1_in + p$h_in)
    w_raw <- p$w_in * 72
    h_raw <- p$h_in * 72
    inner <- sub("^.*?(<g )", "\\1", p$svg_txt)
    inner <- sub("</svg>\\s*$", "", inner)
    sprintf('<svg x="%f" y="%f" width="%f" height="%f" viewBox="0 0 %f %f">\n%s\n</svg>',
            x1_raw, y_top_raw, w_raw, h_raw, p$gv_w_pt, p$gv_h_pt, inner)
  }, character(1))

  svg_txt <- paste(lines_all, collapse = "\n")
  svg_txt <- sub("</svg>\\s*$", paste0(paste(blocks, collapse = "\n"), "\n</svg>\n"), svg_txt)
  writeLines(svg_txt, output, useBytes = TRUE)
  invisible(NULL)
}

# ---- 图例 ----

# 每个type一条竖排dispersion渐变+两端数值(2026-09-13新增)，默认(未命中
# highlight表)桶固定排最前面用default_label当名字(可以是point_label，也
# 可以是"その他"这种更贴切的说法)，之后按highlight表里type出现顺序追加。
.draw_legend <- function(type_ramps, signed_range, legend_cols, default_label, pin,
                          in_to_data_x, in_to_data_y, default_palette = default_color_palette()) {
  legend_w_in <- 0.11
  legend_h_in <- 0.7
  # bar靠右(右边缘固定在离画布右边缘0.05in处，跟下面type标签同一个右对齐
  # 锚点)，数值标注挪到bar左侧、右对齐着往左长；标题也改右对齐，不再居中
  # (2026-09-14确认"图示中的数字应该放在左侧，bar居右，legend的title也
  # 应该居右")
  legend_x2_in <- pin[1] - 0.05
  legend_x1_in <- legend_x2_in - legend_w_in
  # 往下挪，给上方新增的点大小图例(.draw_size_legend)、以及更上方可能有的
  # 形状图例(.draw_shape_legend，合并多种非pair field时才出现)留出空间，
  # 不要贴在一起——1.0改成1.3是2026-09-23加入形状图例后腾出来的额外间距
  # (原1.0在有形状图例时会跟它的最后一项文字撞在一起，实测确认)。
  legend_y1_in <- pin[2] - legend_h_in - 1.3
  # legend_cols是.compute_point_colors算好的现成颜色(线性映射，见
  # .linear_ramp)，这里只管按顺序画矩形条，不重新决定怎么归一化——两端
  # 数字是signed_range里的真实最小/最大值(可以是负的)，颜色沿竖条的
  # 位置也跟这两个数字线性对应(2026-10-17用户确认"直接用归一化"，取值
  # 范围两端数字跟之前一样显示，只是不再强制从0开始)。
  n_steps <- length(legend_cols)
  step_h_in <- legend_h_in / n_steps
  for (i in seq_len(n_steps)) {
    rect(in_to_data_x(legend_x1_in), in_to_data_y(legend_y1_in + (i - 1) * step_h_in),
         in_to_data_x(legend_x2_in), in_to_data_y(legend_y1_in + i * step_h_in),
         col = legend_cols[i], border = NA, xpd = NA)
  }
  rect(in_to_data_x(legend_x1_in), in_to_data_y(legend_y1_in),
       in_to_data_x(legend_x2_in), in_to_data_y(legend_y1_in + legend_h_in),
       border = "#888888", xpd = NA)
  legend_cex <- 0.5
  text(in_to_data_x(legend_x2_in), in_to_data_y(legend_y1_in + legend_h_in),
       "分散度差", cex = legend_cex, adj = c(1, -0.6), xpd = NA)
  text(in_to_data_x(legend_x1_in - 0.03), in_to_data_y(legend_y1_in + legend_h_in),
       format(signed_range[2], digits = 3), cex = legend_cex, adj = c(1, 1), xpd = NA)
  text(in_to_data_x(legend_x1_in - 0.03), in_to_data_y(legend_y1_in),
       format(signed_range[1], digits = 3), cex = legend_cex, adj = c(1, 0), xpd = NA)

  swatch_size_in <- 0.06
  swatch_gap_in <- 0.03
  # type_ramps为空说明没有highlight表/pairs命中，全部点只有default_palette
  # 一种颜色——这种情况下画"その他"色块跟正文没有第二种颜色可区分，纯属
  # 噪音，直接跳过整段色块+文字(2026-09-23确认"如果只有一种颜色的点，
  # legend不需要その他")。
  if (length(type_ramps) == 0) return(invisible(NULL))
  all_labels <- c(list(default_palette), type_ramps)
  names(all_labels)[1] <- default_label
  # type名字长度不固定(比如"近藤(2011)の抽出語の構成要素"比"is_word"长
  # 很多)，从legend_x1_in往右画文字很容易超出画布右边框(2026-09-14发现
  # "legend的文字有些超出了图片的右侧边框")——改成文字右对齐、锚点钉在
  # 离画布右边缘固定一点距离的位置，文字往左长而不是往右长，色块摆在
  # 文字左边(位置跟着文字宽度算，不是固定offset)，这样不管标签多长都
  # 不会超出右边框
  label_right_in <- pin[1] - 0.05
  for (i in seq_along(all_labels)) {
    y_in <- legend_y1_in - swatch_gap_in - (i - 1) * (swatch_size_in + swatch_gap_in) - swatch_size_in
    pal <- all_labels[[i]]
    swatch_col <- colorRampPalette(c(pal["light"], pal["dark"]))(2)[2]
    label <- names(all_labels)[i]
    text(in_to_data_x(label_right_in), in_to_data_y(y_in + swatch_size_in / 2),
         label, cex = legend_cex, adj = c(1, 0.5), xpd = NA)
    label_w_in <- strwidth(label, cex = legend_cex, units = "inches")
    swatch_x2_in <- label_right_in - label_w_in - swatch_gap_in
    swatch_x1_in <- swatch_x2_in - swatch_size_in
    rect(in_to_data_x(swatch_x1_in), in_to_data_y(y_in),
         in_to_data_x(swatch_x2_in), in_to_data_y(y_in + swatch_size_in),
         col = swatch_col, border = NA, xpd = NA)
  }
}

# 点大小代表total_freq(粗频度)，之前一直没有图例说明这层映射
# (2026-09-14确认"还缺少一个，点的大小相关的legend")。挑3个代表性频次
# (min/median/max)，各画一个跟散点图正文同一个size_scale公式算出来的点，
# 数值标在左边(跟dispersion渐变条一样的"数字在左、图形在右"布局)，放在
# 画布右上角、dispersion渐变条上方那一小块预留空间里。
# NODE_SIZEはN_TARGETと同じ順番・同じ長さのベクトル(visprepare.py側で
# 表生成時に書き込み済みのnode_size列、.size_scale自身のコメント参照)。
# 凡例は代表値3つ(最小/中央値/最大のn_target)を見せたいが、中央値は
# 偶数個の時2つの中間値の平均になり、実データのどの行の値とも一致しない
# ことがある -- 「一番近い実際の行」を探してその行自身のn_target(表示
# する数字)とnode_size(.size_scaleにそのまま渡す)を両方読む、Rでlog2/
# maxを計算し直すことはしない(2026-10-17、ユーザー確認済み：表生成時に
# 1回だけ計算)。
.draw_size_legend <- function(n_target, node_size, pin, in_to_data_x, in_to_data_y, pch = 1) {
  legend_cex <- 0.5
  x_in <- pin[1] - 0.15
  target_values <- c(min(n_target), median(n_target), max(n_target))
  closest_idx <- vapply(target_values, function(v) which.min(abs(n_target - v)), integer(1))
  labels <- n_target[closest_idx]
  sizes <- node_size[closest_idx]
  gap_in <- 0.16
  y_top_in <- pin[2] - 0.22
  text(in_to_data_x(pin[1] - 0.05), in_to_data_y(pin[2] - 0.06),
       "頻度", cex = legend_cex, adj = c(1, 0), xpd = NA)
  for (i in seq_along(labels)) {
    y_in <- y_top_in - (i - 1) * gap_in
    points(in_to_data_x(x_in), in_to_data_y(y_in), pch = pch,
           cex = .size_scale(sizes[i]), col = "#666666", xpd = NA)
    text(in_to_data_x(x_in - 0.12), in_to_data_y(y_in), labels[i],
         cex = legend_cex, adj = c(1, 0.5), xpd = NA)
  }
}

# d$point_pch现在可能逐点不同(合并了多种非pair field时)，term frequency
# 这个大小图例本身画的是"形状不重要、只看大小"的3个代表点——形状不统一时
# 挑一个中性的空心圆(pch=1)代表即可，不需要为这一个纯粹讲大小的图例再引入
# 一套"用哪个形状"的判断逻辑。
.representative_pch <- function(pch_vec) {
  u <- unique(pch_vec)
  if (length(u) == 1) u else 1
}

# entry/concept/topic合并展示时(d$field逐行不同、pch也跟着逐点不同)，
# 用一个小图例说明"哪个形状对应哪个field"——只有一种field(没有可比对象)
# 时形状差异没有意义，直接跳过，不画(2026-09-23比照颜色图例"只有一种颜色
# 不需要その他"同一个原则)。放在点大小图例下方、dispersion渐变条上方那
# 一小块预留空间里，同一个x_in锚点对齐、"文字在左、图形在右"同一套布局。
.draw_shape_legend <- function(d, pin, in_to_data_x, in_to_data_y,
                                shapes = field_shape_by_type(), labels = field_label_by_type()) {
  if (is.null(d$field) || !is.null(d$w1)) return(invisible(NULL))
  fields <- sort(unique(d$field))
  if (length(fields) <= 1) return(invisible(NULL))
  legend_cex <- 0.5
  x_in <- pin[1] - 0.15
  gap_in <- 0.14
  y_top_in <- pin[2] - 0.68
  for (i in seq_along(fields)) {
    y_in <- y_top_in - (i - 1) * gap_in
    pch_i <- unname(shapes[fields[i]])
    if (is.na(pch_i)) pch_i <- 16
    points(in_to_data_x(x_in), in_to_data_y(y_in), pch = pch_i, cex = 1, col = "#666666", xpd = NA)
    label_i <- unname(labels[fields[i]])
    if (is.na(label_i)) label_i <- fields[i]
    text(in_to_data_x(x_in - 0.12), in_to_data_y(y_in), label_i,
         cex = legend_cex, adj = c(1, 0.5), xpd = NA)
  }
}

# Replaces an earlier .draw_end_labels (removed 2026-10-16) that drew
# value_label/not_value_label as large text INSIDE the data panel's
# own corners, which collided with points/badges. This instead draws
# one real double-headed arrow (arrows(code=3), not a plotmath %<->%
# character -- user explicitly wants a drawn line, not a glyph, 2026-
# 10-16) per axis, spanning that axis's own FULL length (usr[1] to
# usr[2] for X, usr[3] to usr[4] for Y -- "宽度和轴的长度相同"), placed
# in the margin between the tick labels and that axis's own title.
# Margin "line" units are converted to data coordinates via par("csi")
# (inches per line, independent of cex/canvas_scale) and the panel's
# own inches-per-data-unit ratio (same ratio panel_coord_converters()
# computes, inlined here since this runs in the margin, outside usr,
# where panel_coord_converters's own in_to_data_*(positive inches)
# convention isn't what's needed -- this wants negative/leftward
# offsets instead).
# X axis keeps value_label/not_value_label (the point->target KLD's
# own sign, unchanged meaning). Y axis is "uniform"<->"concentrated" --
# dispersion's own meaning (low = concentrated in a few texts =
# "concentrated", high = spread evenly relative to those texts' own
# weight = "uniform"), fixed wording, not a caller-supplied label,
# since dispersion's semantics don't vary by dataset the way
# value_label does -- "concentrated" goes at the bottom (low
# dispersion, near the axis origin), "uniform" at the top (high
# dispersion). Renamed 2026-10-17 from "specific"/"common": dispersion
# only measures how proportionally a point's occurrences are spread
# across the handful of lines it itself touches (see
# level2_line_matrix's own docstring -- its population total is scoped
# to those lines alone, never the whole corpus), so it carries no
# information about how widely/frequently a word is actually used.
# "common" wrongly implied the latter, and is actively misleading for
# low-frequency points: a hapax (n=1 line touched) always gets
# dispersion EXACTLY 1.0 by construction (P and Q collapse to the same
# single-point distribution, giving KLD=0), regardless of whether the
# word is remotely common in any ordinary sense -- confirmed with the
# user via a concrete example (あかしがた【明石潟】, total_freq=1,
# dispersion=1.0). "uniform"/"concentrated" describe the actual
# mechanism (proportional-to-local-weight vs skewed) without that
# false implication. The Y axis itself is plotted the ordinary way
# (usr[3]=min/bottom, usr[4]=max/top -- the 2026-10-15 axis swap moved
# dispersion onto Y but didn't flip its own orientation), so this has
# to match _normalize_dispersion's own direction (close to 0 =
# concentrated, close to 1 = uniform -- see that function's own
# docstring in visprepare.py) directly, not the reverse.
.draw_axis_direction_arrows <- function(value_label, not_value_label) {
  usr <- par("usr"); pin <- par("pin")
  x_per_in <- diff(usr[1:2]) / pin[1]
  y_per_in <- diff(usr[3:4]) / pin[2]
  line_in <- par("csi")
  mgp <- par("mgp")
  # Three separate rows now, each its own margin line, outward from the
  # plot: tick labels (occupy roughly [mgp[2], mgp[2]+1]) -> arrow (at
  # arrow_line) -> end-text labels (at label_line, a further 0.5 lines
  # out so the text doesn't sit on top of the arrowhead) -> axis title
  # (furthest out, at .draw_axes's own title_line = mgp[1] + 1.4) --
  # see .open_canvas's own comment for the matching mar increase that
  # makes room for all three new rows.
  arrow_line <- mgp[2] + 1.4
  label_line <- arrow_line + 0.5
  label_cex <- 0.8

  x_arrow_y <- usr[3] - arrow_line * line_in * y_per_in
  x_label_y <- usr[3] - label_line * line_in * y_per_in
  arrows(usr[1], x_arrow_y, usr[2], x_arrow_y, code = 3, length = 0.05,
         angle = 25, col = "#888888", lwd = 0.8, xpd = NA)
  text(usr[1], x_label_y, not_value_label, adj = c(0, 0.5), cex = label_cex,
       col = "#888888", xpd = NA)
  text(usr[2], x_label_y, value_label, adj = c(1, 0.5), cex = label_cex,
       col = "#888888", xpd = NA)

  y_arrow_x <- usr[1] - arrow_line * line_in * x_per_in
  y_label_x <- usr[1] - label_line * line_in * x_per_in
  arrows(y_arrow_x, usr[3], y_arrow_x, usr[4], code = 3, length = 0.05,
         angle = 25, col = "#888888", lwd = 0.8, xpd = NA)
  text(y_label_x, usr[3], "concentrated", srt = 90, adj = c(0, 0.5), cex = label_cex,
       col = "#888888", xpd = NA)
  text(y_label_x, usr[4], "uniform", srt = 90, adj = c(1, 0.5), cex = label_cex,
       col = "#888888", xpd = NA)
}

# ---- draw_point_vs_target()自己的准备步骤 ----

resolve_annotation_inputs <- function(annotate_words, annotate_words_file,
                                       highlight, highlight_file,
                                       highlight_pairs, highlight_pairs_file,
                                       annotate_pairs = NULL, annotate_pairs_file = NULL) {
  if (!is.null(annotate_words_file)) {
    raw <- readLines(annotate_words_file, encoding = "UTF-8")
    annotate_words <- unique(raw[nchar(raw) > 0])
  }
  if (!is.null(highlight_file)) {
    highlight <- read.table(highlight_file, sep = "\t", header = TRUE,
                             stringsAsFactors = FALSE, quote = "", fileEncoding = "UTF-8")
  }
  if (!is.null(highlight_pairs_file)) {
    highlight_pairs <- read.table(highlight_pairs_file, sep = "\t", header = TRUE,
                                    stringsAsFactors = FALSE, quote = "", fileEncoding = "UTF-8")
  }
  if (!is.null(annotate_pairs_file)) {
    annotate_pairs <- read.table(annotate_pairs_file, sep = "\t", header = TRUE,
                                  stringsAsFactors = FALSE, quote = "", fileEncoding = "UTF-8")
  }
  list(annotate_words = annotate_words, highlight = highlight, highlight_pairs = highlight_pairs,
       annotate_pairs = annotate_pairs)
}

# Named by every column either the word- or pair-level visualization
# tsv schema can carry (see load_and_merge_input's own fread() call
# below) -- fread only applies colClasses entries whose name matches
# an actual column in THIS file; names from the OTHER schema just
# aren't present in any one file, which fread warns about (benign,
# suppressWarnings()'d below) rather than erroring on.
.POINT_TABLE_COL_CLASSES <- c(
  point = "character", frequency_type = "character", kld_variant = "character",
  keyness = "numeric", frequency = "numeric", total_freq = "numeric",
  n_target = "numeric", n_not_target = "numeric", node_size = "numeric",
  dispersion = "numeric", dispersion_difference = "numeric", field = "character",
  w1 = "character", w2 = "character", field1 = "character", field2 = "character",
  freq1 = "numeric", freq2 = "numeric", dispersion1 = "numeric", dispersion2 = "numeric",
  assoc_w1_to_w2 = "numeric", assoc_w2_to_w1 = "numeric"
)

load_and_merge_input <- function(input) {
  # inputはカンマ区切りで複数パス指定可——pair-entry-entry/pair-entry-concept/
  # pair-entry-topicのような複数のenriched tsvを1枚の図にまとめて描く時に使う。
  # 列が完全に揃っていない場合(field1/field2がまだ無い古いキャッシュtsvなど)
  # でも落ちないよう.rbind_fill()で欠けている列はNA埋めする——単一ファイル
  # 指定という既存の使い方は今までどおり(splitしても要素数1、rbind_fillも
  # 1個だけならそのまま返す)。
  #
  # data.table::fread(), not read.table() -- 2026-10-17, confirmed
  # against the real pair-entry-entry-vs-anthology-shinkokin.tsv
  # (1.8M rows): read.table (even tuned with encoding=/colClasses) was
  # still ~5s on that file, too slow for Shiny's interactive adjust-
  # and-preview loop; fread is built for exactly this. data.table =
  # FALSE -- every downstream filter_by_*/subset in this file does
  # plain d[keep, ] data.frame-style row filtering, which a
  # data.table's own `[` would interpret with different (NSE)
  # semantics; converting back to a plain data.frame right at read
  # time keeps every existing call site correct unchanged.
  input_paths <- strsplit(input, ",", fixed = TRUE)[[1]]
  tables <- lapply(input_paths, function(p) {
    # Keep existing TSV commands usable with the compact published assets.
    rds_path <- if (grepl("\\.rds$", p)) p else sub("\\.tsv$", ".rds", p)
    if (!file.exists(p) && !file.exists(rds_path) &&
        exists(".base_url", envir = .GlobalEnv, inherits = FALSE)) {
      .fetch_browser_asset(paste0("data/processed/visualization/", basename(rds_path)), rds_path)
    }
    if (grepl("\\.rds$", p) || (!file.exists(p) && file.exists(rds_path))) {
      d <- readRDS(rds_path)
      stopifnot(is.data.frame(d), "point" %in% names(d))
      return(d)
    }
    suppressWarnings(data.table::fread(
      p, sep = "\t", header = TRUE, quote = "", encoding = "UTF-8",
      na.strings = c("NA", ""), colClasses = .POINT_TABLE_COL_CLASSES,
      data.table = FALSE
    ))
  })
  if (length(tables) == 1) tables[[1]] else .rbind_fill(tables)
}

resolve_zoom_box <- function(zoom_box, zoom_x1, zoom_x2, zoom_y1, zoom_y2, badge_selection, d) {
  # zoom_x1/x2/y1/y2是zoom_box的CLI友好写法(4个边界各自可选，没给的边界
  # 退回数据本身的min/max)。badge_selection="zoom_box"时即使4个边界都没给
  # 也必须建出zoom_box(全部退回min/max，等价于"整个图的范围都算")，不然
  # .select_badge_rows()的zoom_box分支会因stopifnot(!is.null(zoom_box))报错
  needs_box <- is.null(zoom_box) && (identical(badge_selection, "zoom_box") ||
                                      !is.null(zoom_x1) || !is.null(zoom_x2) ||
                                      !is.null(zoom_y1) || !is.null(zoom_y2))
  if (!needs_box) return(zoom_box)
  list(
    x1 = if (!is.null(zoom_x1)) as.numeric(zoom_x1) else min(d$frequency),
    x2 = if (!is.null(zoom_x2)) as.numeric(zoom_x2) else max(d$frequency),
    y1 = if (!is.null(zoom_y1)) as.numeric(zoom_y1) else min(d$keyness),
    y2 = if (!is.null(zoom_y2)) as.numeric(zoom_y2) else max(d$keyness)
  )
}

# resolve_zoom_box()算出来的边界是pre-break坐标(用户在轴上读到/CLI填的
# 真实KLD值，或者pre-break d的min/max兜底)，但断轴一旦触发，从这之后
# d$frequency/d$keyness全部换成了压缩后的坐标——zoom_box不跟着换算的话，
# 就是在拿"真实KLD值"当"压缩后坐标"直接用，两者数量级往往差几十倍
# (fwd()在主区间的斜率是main_frac/main_end，通常远大于1)，画出来的框会
# 紧贴轴根部退化成一条线，badge_selection="zoom_box"的筛选也会连带出错
# (2026-09-25用户反馈"扩大表示时主图里的框横轴范围变得很奇怪，没有
# 宽度")。fwd()是单调的，min/x1/y1和max/x2/y2各自单独换算，结果仍然
# 是x1<x2、y1<y2的合法矩形，不需要额外排序。zoom_box为NULL(没用zoom_box
# 相关选项)时原样返回，不受影响。
.transform_zoom_box_for_breaks <- function(zoom_box, breaks) {
  if (is.null(zoom_box)) return(zoom_box)
  if (!is.null(breaks$x_break)) {
    zoom_box$x1 <- breaks$x_break$fwd(zoom_box$x1)
    zoom_box$x2 <- breaks$x_break$fwd(zoom_box$x2)
  }
  if (!is.null(breaks$y_break)) {
    zoom_box$y1 <- breaks$y_break$fwd(zoom_box$y1)
    zoom_box$y2 <- breaks$y_break$fwd(zoom_box$y2)
  }
  zoom_box
}

compute_point_shapes <- function(d, point_field) {
  # is_pair版(有w1列)维持"整张图统一一种形状"(pair合并展示的形状区分只放
  # 在network小图里)；非pair(word)版合并了多种field(entry/concept/topic)时
  # (d$field逐行不同)逐点套用不同形状；没有d$field(老缓存/未合并场景)退回
  # 原来"entry圆、非entry方"整图统一一种形状
  if (is.null(d$w1) && !is.null(d$field)) {
    pch_vec <- unname(field_shape_by_type()[d$field])
    pch_vec[is.na(pch_vec)] <- 1
    pch_vec
  } else {
    rep(if (identical(point_field, "entry")) 1 else 0, nrow(d))
  }
}

select_badges <- function(d, annotate_words, label_top_n, neutral, annotate_nontarget_only,
                           annotate_point_exclude, badge_selection, frontier_n, zoom_box,
                           aggregate_locators, annotate_pairs = NULL) {
  rows <- .select_badge_rows(d, annotate_words, label_top_n, neutral,
                              annotate_nontarget_only, annotate_point_exclude,
                              badge_selection = badge_selection, frontier_n = frontier_n,
                              zoom_box = zoom_box, annotate_pairs = annotate_pairs)
  # 同一locator(n_target/n_not_target完全相同、同一种pair/condition)上堆了
  # 不止1个候选时合并成一个badge，见.aggregate_locator_groups()自己的注释；
  # aggregate_locators=FALSE时跳过，退回聚合功能加入之前的行为
  if (aggregate_locators) rows <- .aggregate_locator_groups(rows)
  rows
}

# usr/pin(当前绘图区的坐标范围+物理英寸尺寸)以及跟着它俩算出来的
# data<->inches双向换算闭包——非inset主图、inset上下两块面板，一共4处
# 分别调用，每次都是同一套算法(读当前par()、算换算系数)，只是各自面板的
# usr/pin不同。全程在inches空间做badge排版/连线折角是因为x/y轴数量级往往
# 差好几个量级，直接用data坐标算距离/角度会被轴尺度扭曲。
panel_coord_converters <- function() {
  usr <- par("usr")
  pin <- par("pin")
  x_per_in <- (usr[2] - usr[1]) / pin[1]
  y_per_in <- (usr[4] - usr[3]) / pin[2]
  list(
    usr = usr, pin = pin, x_per_in = x_per_in, y_per_in = y_per_in,
    data_to_in_x = function(x) (x - usr[1]) / x_per_in,
    data_to_in_y = function(y) (y - usr[3]) / y_per_in,
    in_to_data_x = function(x) usr[1] + x * x_per_in,
    in_to_data_y = function(y) usr[3] + y * y_per_in
  )
}

# ---- 对外入口 ----

# 两份独立的名单(2026-09-13确认"命名成两种列表"，取代之前合在一起的
# annotations表)，可以直接传内存对象，也可以传文件路径现读(2026-09-13
# 确认"直接写到文件，增加两个读取关联文件的参数"，见
# scripts/analysis/prep_kondo_lists.R这个一次性预处理脚本)：
#   annotate_words/annotate_words_file -- 单列的字符向量/单列文本文件
#     (一行一个読み字符串)，命中的点画完整badge，不带类别；badge本身的
#     文字已经是它的视觉信号，不需要类别配色。
#   highlight/highlight_file -- 两列的data.frame/TSV(term, type)，命中的
#     点按type分类别染色(不画badge)，每个type领一条独立渐变色带、进图例。
#   两个_file参数如果给了，会覆盖对应的直接传值参数。
# annotate_nontarget_only/annotate_point_exclude只作用在annotate_words筛出来的
# 候选上，分别对应"只保留某种收窄条件(比如n_target==0)的点"和"手动排除
# 同音异字误伤"。
# default_label：没被highlight命中的点在图例里用什么名字(默认沿用
# point_label，比如"その他"这种更贴切的说法可以自己传)。
# badge_layout: "empty_region"(新默认，找最大空白矩形)或"repel"(力导向，
# 重新引入ggrepel)。
# zoom_annotated_only/repel_annotated_only：两个彼此独立的"只画被标注的点"
# 开关(2026-09-24"扩大局部时可以选择不展示不被annotate的点"；后续追加
# "repel模式下也要有类似的参数，但要区分成两个函数")。前者只管inset的
# 扩大面板，后者只要badge_layout=="repel"就生效、不分主图还是扩大面板——
# 两者是"或"的关系，扩大面板里两个条件命中任一个都会精简；共用同一个
# filter_to_annotated_points()筛选原语，判断逻辑本身分成
# should_hide_unannotated_in_zoom()/should_hide_unannotated_for_repel()
# 两个独立函数，不是同一个参数复用。默认都是FALSE，不开启时行为不变。
# inset: 开关"扩大图"(2026-09-15确认最终设计，取代之前usage.txt里嵌入式
# 子图的草案，以及本次会话中途做过的两版中间实现——依次被用户否掉，理由
# 分别是"应该是独立的完整图片，不是嵌在主图里面"、"badge应该画在扩大图里，
# 不是主图"、"框体没有对齐四个角和连线没有结合"+"框需要有margin"(这一版
# 一度把框做成"比取景范围小一圈"，后来发现理解错了——见下面inset_pad_frac
# 的说明)。默认FALSE，行为完全不变(见下面if(inset)分支，非inset路径原样
# 保留，一行没动)。inset=TRUE时画布整体加宽一倍(canvas_size不变，宽度变成
# 2倍)，左边是主图(不画badge，只画一个虚线框标出被放大的区域)，右边是
# 一张独立、带完整坐标轴的正方形放大图(badge画在这里)，右图额外用box()
# 画一圈四条边都有的实线红框(主图惯用的"只画左/下轴线"在这里不够，右图
# 要看起来像个完整独立的画框)，两图之间用虚线(跟已有neutral/median参考线
# 的线型区分开)从框的四个角精确连到框的四个角，仿照参考PDF(Nature复现
# 文章)里ggmagnify的"projection"效果。
# 放大范围不是手动指定的xlim/ylim，而是自动取badge候选点(annotate_words/
# frontier选中的那批词，不是全部点)的外接框("扩大的范围是annotate的词所
# 被囊括的...区块")，这个框同时是(a)主图里标出来的那个框、也是(b)放大图
# 自己的坐标轴显示范围，两边共用同一个矩形，框本身就是取景范围，不再
# 额外留白。inset_pad_frac是这个框往外扩的留白比例——它解决的是"badge
# 候选点在某个轴上数值高度集中(比如kokin最上面那批词keyness几乎相同)时
# 外接框在那个维度上会退化成一条线"，不是"框跟坐标轴之间留白"(2026-09-15
# 确认"margin的意思是保证退化维度有个像样的最小尺寸，不是框比取景范围
# 小"这点反馈后修正)。inset_border_col是两边的框+连线的颜色。

# 扩大面板(render_zoomed_scatter_panel)专用：ZOOM_ANNOTATED_ONLY为TRUE时，
# 扩大框范围内没被选进badge_rows的点不画，跟badge_layout用哪种布局无关。
should_hide_unannotated_in_zoom <- function(zoom_annotated_only) {
  isTRUE(zoom_annotated_only)
}

# 主图/扩大面板都适用：REPEL_ANNOTATED_ONLY为TRUE且badge_layout=="repel"时，
# 未被标注的点不画——只在repel布局下有意义(empty_region布局找空白矩形，
# 不需要"点越少越好摆"这个动机)，因此额外检查badge_layout，不是
# repel_annotated_only参数本身就够。
should_hide_unannotated_for_repel <- function(repel_annotated_only, badge_layout) {
  isTRUE(repel_annotated_only) && identical(badge_layout, "repel")
}

draw_point_vs_target <- function(input, output, point_label = "point",
                                  axis_label = "target", value_label = "value",
                                  not_value_label = paste0("¬", value_label),
                                  label_top_n = 10, nofreq = FALSE,
                                  annotate_words = NULL,
                                  annotate_words_file = NULL,
                                  highlight = NULL,
                                  highlight_file = NULL,
                                  highlight_pairs = NULL,
                                  highlight_pairs_file = NULL,
                                  annotate_nontarget_only = FALSE,
                                  annotate_point_exclude = NULL,
                                  highlight_frontier = FALSE,
                                  default_label = point_label,
                                  badge_layout = "empty_region",
                                  badge_selection = "words",
                                  frontier_n = 15,
                                  canvas_size = 450,
                                  inset = FALSE,
                                  inset_position = "top",
                                  inset_border_col = "#FF073A",
                                  inset_pad_frac = 0.15,
                                  zoom_annotated_only = FALSE,
                                  repel_annotated_only = FALSE,
                                  zoom_box = NULL,
                                  zoom_x1 = NULL,
                                  zoom_x2 = NULL,
                                  zoom_y1 = NULL,
                                  zoom_y2 = NULL,
                                  edge_threshold = NULL,
                                  dispersion_threshold = NULL,
                                  node_dispersion_threshold = NULL,
                                  frequency_threshold = NULL,
                                  node_frequency_threshold = NULL,
                                  aggregate_locators = TRUE,
                                  max_group_members = 50,
                                  network_scale = 0.01,
                                  network_node_gap = 1,
                                  point_field = "entry",
                                  annotate_pairs = NULL,
                                  annotate_pairs_file = NULL,
                                  kld_variant = NULL,
                                  frequency_type = NULL) {
  # badge_layout="network"splices graphviz output in as svg XML text
  # (.splice_network_glyphs below) -- that only makes sense for a real
  # vector output, never png. Rather than make the caller get this
  # right and stop() otherwise, just force the extension here, before
  # OUTPUT is used anywhere below (both the inset and non-inset paths
  # read this same local variable).
  if (identical(badge_layout, "network") && !grepl("\\.svg$", output, ignore.case = TRUE)) {
    output <- sub("\\.[^.]*$", ".svg", output)
  }
  inputs <- resolve_annotation_inputs(annotate_words, annotate_words_file,
                                       highlight, highlight_file,
                                       highlight_pairs, highlight_pairs_file,
                                       annotate_pairs, annotate_pairs_file)
  annotate_words <- inputs$annotate_words
  highlight <- inputs$highlight
  highlight_pairs <- inputs$highlight_pairs
  annotate_pairs <- inputs$annotate_pairs
  d <- load_and_merge_input(input)
  # scripts/patches/visprepare.py now always merges all of its own
  # POINT_VARIANTS ("obs"/"nofreq"/"norm") into one tsv with a
  # "kld_variant" column rather than writing one tsv per variant --
  # KLD_VARIANT, if given, is this function's own single place to pick
  # one variant's rows back out (trust upstream shape: the column is
  # assumed present whenever this is non-NULL; NULL is a complete
  # no-op, so any older tsv without the column -- or a non-point
  # reference tsv -- is unaffected).
  if (!is.null(kld_variant)) d <- d[d$kld_variant == kld_variant, ]
  # frequency_type works the same way, one axis over -- visprepare.py's
  # own FREQUENCY_TYPES ("term"/"doc") are merged into the SAME tsv
  # too, unlike kld_variant this axis also changes n_target/n_not_
  # target/total_freq (and, for pairs, assoc_w1_to_w2/assoc_w2_to_w1),
  # not just keyness/frequency -- see visprepare.py's own module
  # docstring. Still just one more row filter here either way.
  if (!is.null(frequency_type)) d <- d[d$frequency_type == frequency_type, ]
  edge_threshold <- if (!is.null(edge_threshold)) as.numeric(edge_threshold) else NULL
  dispersion_threshold <- if (!is.null(dispersion_threshold)) as.numeric(dispersion_threshold) else NULL
  node_dispersion_threshold <- if (!is.null(node_dispersion_threshold)) as.numeric(node_dispersion_threshold) else NULL
  frequency_threshold <- if (!is.null(frequency_threshold)) as.numeric(frequency_threshold) else NULL
  node_frequency_threshold <- if (!is.null(node_frequency_threshold)) as.numeric(node_frequency_threshold) else NULL
  d <- filter_by_edge_strength(d, edge_threshold)
  d <- filter_by_dispersion(d, dispersion_threshold)
  d <- filter_by_node_dispersion(d, node_dispersion_threshold)
  d <- filter_by_frequency(d, frequency_threshold)
  d <- filter_by_node_frequency(d, node_frequency_threshold)
  # Thresholds set too strict can filter every row out, which later makes
  # range(d$frequency)/range(d$keyness) return c(Inf, -Inf) and crash deep
  # inside plot.window() with a cryptic "need finite 'xlim' values" -- fail
  # here instead, right where the actual cause (an empty d) is known.
  if (nrow(d) == 0) {
    stop("No points left after filtering -- edge_threshold/dispersion_threshold/",
         "node_dispersion_threshold/frequency_threshold/node_frequency_threshold ",
         "are too strict for this dataset")
  }
  # Snapshot taken before the frequency/keyness overwrite below -- the
  # annotated-points table (app.R's table tab) wants the real tsv-shaped
  # columns, not the plotting-convenience renames. .row_id survives
  # rbind()/subsetting downstream even when row names don't (rbind()
  # silently renames colliding ones via make.unique), so it is the robust
  # join key back from tie_box_rows to this snapshot, done once at the
  # very end of this function.
  d$.row_id <- seq_len(nrow(d))
  d_raw_cols <- d
  # 纵軸→dispersion、横軸→point2target(旧keyness)に恒久的に差し替え
  # (2026-10-15、ユーザー確認済み：frequencyは図から完全に外す。トグルで
  # はなく今後はこれが唯一のビュー)。全てのフィルタ(上のfilter_by_*)は
  # まだ元のdispersion/total_freq/keyness列を読む必要があるので、ここ
  # (フィルタの後、zoom_box/ranges計算の前)が差し替えの場所——下流の
  # zoom_box解決/badge選定/ネットワークグリフ配置は全部d$frequency(x)/
  # d$keyness(y)を素直に読むだけなので、ここで値を入れ替えておけば無改修
  # で動く。neutral参考線/hi-lo色分け(どちらもd$keyness>neutralという
  # 符号付きKLD前提のロジック)はdispersionには意味を持たない組み合わせに
  # なるが、ユーザー確認済み：今は直さず既知の見た目のクセとして残す。
  d$frequency <- d$keyness
  d$keyness <- d$dispersion
  # zoom_box在这里(用pre-break的d$frequency/keyness当min/max兜底)解析，
  # 早于下面的.apply_nofreq_breaks——2026-09-22发现的bug("选了zoom_box选项
  # 但还没填4个边界时会报错")修复后就是这个顺序，沿用不变
  zoom_box <- resolve_zoom_box(zoom_box, zoom_x1, zoom_x2, zoom_y1, zoom_y2, badge_selection, d)

  ranges <- .compute_axis_ranges(d, nofreq)
  breaks <- .apply_nofreq_breaks(d, nofreq)
  d <- breaks$d
  # badge筛选(hi/lo按keyness跟neutral比较)要用跟d$keyness同一个坐标系——
  # d$keyness如果被断轴压缩过，neutral(理论中性点)也要用同一个fwd()换算
  neutral <- if (!is.null(breaks$y_break)) breaks$y_break$fwd(ranges$neutral) else ranges$neutral
  # zoom_box跟上面的neutral同一类问题(2026-09-25用户反馈"扩大表示时主图里
  # 的框横轴范围变得很奇怪，没有宽度")：zoom_box是用pre-break的d算出来的
  # (见上面resolve_zoom_box那行的注释)，取的是用户在图上读到的真实KLD值，
  # 但从这一行往下badge_rows/.draw_source_box/render_zoomed_scatter_panel
  # 用的d$frequency/d$keyness都已经是断轴压缩后的坐标——不换算的话，
  # zoom_box的边界值(比如0.001~0.003这种远小于压缩后典型取值的原始KLD数)
  # 会被当成压缩后坐标直接使用，画出来的框自然退化成贴在轴根部的一条线，
  # badge_selection="zoom_box"的取景筛选也会因此选中几乎全部或几乎没有点。
  zoom_box <- .transform_zoom_box_for_breaks(zoom_box, breaks)
  # 颜色要在选badge_rows之前算好、点到d$point_col上——badge_rows是d的一份
  # 子集快照(data.frame子集是复制，不是引用)，先选badge_rows、后面才在d上
  # 补point_col列的话，badge_rows拿到的还是没有point_col列的旧快照
  # (2026-09-15曾经踩过的bug："badge在附图里"看起来不对，根源在此)
  colors <- .compute_point_colors(d, highlight, highlight_pairs, annotate_point_exclude,
                                   highlight_frontier, frontier_n)
  d$point_col <- colors$point_col
  # 存成d的一列而不是独立变量，让它跟着d一起被badge_rows/dd这些子集操作
  # 自动带过去，不用在每处子集代码另外重新计算
  d$point_pch <- compute_point_shapes(d, point_field)
  # badge_rows的选取逻辑不管inset开不开都一样，只是inset=TRUE时它还多了
  # 一个用途(决定放大范围)——提到if/else分支之前算一次，两边共用
  badge_rows <- select_badges(d, annotate_words, label_top_n, neutral,
                               annotate_nontarget_only, annotate_point_exclude,
                               badge_selection, frontier_n, zoom_box, aggregate_locators,
                               annotate_pairs = annotate_pairs)
  # .draw_x_tie_boxes自身のdocstring参照：aggregate_locators=TRUEの時、
  # badge_rowsは既に(n_target, n_not_target, kind)で1行にまとめられた後
  # (.aggregate_locator_groups)なので、同じキーでもう一度グループ化しよ
  # うとしても各グループが必ず1行になり、箱が一つも描けなくなる
  # (2026-10-17にユーザー報告で発覚)。aggregate_locators=FALSEで同じ
  # select_badgesをもう一度呼び、集約前の個々の点(「annotateされた点」
  # ではあるが、まだ合算されていない実際の座標を持つ行)を別に用意する。
  tie_box_rows <- select_badges(d, annotate_words, label_top_n, neutral,
                                 annotate_nontarget_only, annotate_point_exclude,
                                 badge_selection, frontier_n, zoom_box,
                                 aggregate_locators = FALSE,
                                 annotate_pairs = annotate_pairs)
  # The table-tab's "currently annotated points" output -- same row set as
  # tie_box_rows, joined back to the pre-overwrite snapshot (see .row_id's
  # own comment above) so callers see the real tsv columns.
  annotated_table <- d_raw_cols[d_raw_cols$.row_id %in% tie_box_rows$.row_id, , drop = FALSE]
  annotated_table$.row_id <- NULL
  annotated_table <- .format_display_columns(annotated_table)

  if (inset) {
    .render_with_inset(d, badge_rows, tie_box_rows, colors, output, ranges, breaks, neutral, nofreq,
                        axis_label, point_label, value_label, not_value_label,
                        default_label, badge_selection, badge_layout,
                        canvas_size, inset_border_col, inset_pad_frac,
                        inset_position = inset_position,
                        zoom_box = zoom_box, max_group_members = max_group_members,
                        network_scale = network_scale,
                        network_node_gap = network_node_gap,
                        zoom_annotated_only = zoom_annotated_only,
                        repel_annotated_only = repel_annotated_only,
                        kld_variant = kld_variant)
    cat(output, ":", nrow(d), point_label, "points\n")
    return(invisible(annotated_table))
  }

  .open_canvas(output, canvas_size)
  on.exit(dev.off())
  .draw_axes(d, ranges, nofreq, axis_label, point_label, breaks$x_break, breaks$y_break,
             kld_variant = kld_variant)

  # d_drawn只影响画哪些点/badge布局要避开哪些点，图例(.draw_size_legend/
  # .draw_shape_legend)仍然用完整的d——图例描述的是整份数据的编码规则，
  # 不是"这张图实际画了哪些点"
  d_drawn <- if (should_hide_unannotated_for_repel(repel_annotated_only, badge_layout)) {
    filter_to_annotated_points(d, badge_rows)
  } else d
  points(d_drawn$frequency, d_drawn$keyness, pch = d_drawn$point_pch,
         cex = .size_scale(d_drawn$node_size), col = d_drawn$point_col)

  pc <- panel_coord_converters()
  badge_cex <- 0.55
  fam <- getOption("kld.plot.font_family", "Hiragino Mincho ProN")
  label_gap_in <- strheight("あ", cex = badge_cex, family = fam) * 0.5

  tie_boxes <- .compute_x_tie_boxes(
    tie_box_rows, point_radius_in = .cex_to_radius_in(.size_scale(tie_box_rows$node_size))
  )
  .draw_x_tie_boxes(tie_boxes)
  badge_colors <- .badge_color_from_point(badge_rows)
  .draw_badges(badge_rows, fill = badge_colors$fill, border = badge_colors$border_col,
               text_col = badge_colors$text_col, word_col = badge_colors$word_col,
               badge_layout = badge_layout, d_all = d_drawn,
               value_label = value_label, not_value_label = not_value_label,
               badge_cex = badge_cex, fam = fam, label_gap_in = label_gap_in,
               x_per_in = pc$x_per_in, y_per_in = pc$y_per_in,
               data_to_in_x = pc$data_to_in_x, data_to_in_y = pc$data_to_in_y,
               in_to_data_x = pc$in_to_data_x, in_to_data_y = pc$in_to_data_y,
               pin = pc$pin, neutral = neutral, max_group_members = max_group_members,
               tie_boxes = tie_boxes)

  .draw_axis_direction_arrows(value_label, not_value_label)
  .draw_size_legend(d$n_target, d$node_size, pc$pin, pc$in_to_data_x, pc$in_to_data_y,
                     pch = .representative_pch(d$point_pch))
  .draw_shape_legend(d, pc$pin, pc$in_to_data_x, pc$in_to_data_y)
  .draw_legend(colors$type_ramps, colors$signed_range, colors$legend_cols, default_label, pc$pin,
               pc$in_to_data_x, pc$in_to_data_y, colors$default_palette)

  cat(output, ":", nrow(d), point_label, "points\n")
  invisible(annotated_table)
}

# ---- 扩大图(inset) ----
# 2026-09-15确认最终设计，两次否掉之前的方案后定下来的：不是叠在主图内部
# 的小窗口(那样"扩大"的视觉效果不够，坐标轴/badge都施展不开)，而是画布
# 整体加宽一倍——左边主图(不画badge)，右边一张独立、带完整坐标轴的正方形
# 放大图(badge画在这里)，两边之间用虚线框+四角连线做"projection"提示，
# 参考PDF(Nature复现文章)里ggmagnify的效果。

# badge_rows(annotate/frontier选中的那批候选词)的原始数据点(frequency,
# keyness)外接框，按pad_frac往外扩一点留白——这个框同时是(a)主图里标出
# "这块被放大了"的那个框，也是(b)放大图自己的坐标轴显示范围，两边共用
# 同一个矩形，不再分box/view两层(2026-09-15第一版拿"box比view小一圈"
# 当margin，被用户否掉："margin"真正要解决的问题是另一个——badge候选点
# 在某个轴上数值高度集中(比如kokin最上面那批女性词汇keyness几乎相同)时，
# 外接框在那个维度上会退化成一条看不出宽度的线，需要的是"保证退化维度
# 有个像样的最小尺寸"，不是"框比取景范围小")。
# 用pmax(...)而不是"if(diff==0)才退化"：即使两个候选点的y只是挨得很近
# 但不是严格相等(diff(y_rng)是个很小的正数)，算出来的pad仍然会小到肉眼
# 看不出宽度——用主图该轴总range的5%当最小值兜底，不管diff(y_rng)是0
# 还是只是很小，都能保证画出来是个矩形而不是一条线。
.compute_zoom_box <- function(badge_rows, main_ranges, pad_frac, min_span_frac = 0.05) {
  x_rng <- range(badge_rows$frequency); y_rng <- range(badge_rows$keyness)
  xpad <- max(diff(x_rng) * pad_frac, diff(main_ranges$xr) * min_span_frac)
  ypad <- max(diff(y_rng) * pad_frac, diff(main_ranges$yr) * min_span_frac)
  list(x1 = x_rng[1] - xpad, x2 = x_rng[2] + xpad,
       y1 = y_rng[1] - ypad, y2 = y_rng[2] + ypad)
}

# 主图里标出"这块被放大了"的框：虚线(跟已有neutral/median参考线的线型
# 区分开)，只占主图的一小块区域，所以用rect()在data坐标里手画。
.draw_source_box <- function(box, border_col) {
  rect(box$x1, box$y1, box$x2, box$y2, border = border_col, lty = "1343",
       lwd = 1.2, xpd = NA)
}

# 2026-10-17、ユーザー確認済み：X軸(point→target KLD)の値が完全一致する
# tieでグループ化していた最初の版は、低頻度語の大半を取りこぼしていた
# ことが判明した -- pair_cross_table自身のdocstring通り、population
# (分母側)はMEASURED自身が実際に現れる行ごとに計算される(scope全体で
# 共有される定数ではない)ため、n_target/n_not_targetが同じ組み合わせ
# でも実際に現れる行が違えばKLDはわずかに(だが非ゼロに)ズレる--図の
# 上では同じ場所に密集して見えても、浮動小数点としては別の値になる。
# なのでグループ化の基準をX(keyness)の完全一致から、ROWSが実際に持つ
# n_target/n_not_target(target/neg_target側のmatch数、こちらは本当に
# 同じ値になる)の一致に変更。X・Y両方ともその組の値の範囲(min~max)を
# 取った外接矩形を描く(前版は常にX固定の細い縦長の箱だったが、今は同じ
# グループ内でもXが微妙にばらつくため、Xにも幅を持たせる必要がある)。
# MIN_GROUP_SIZE=2(2026-10-17、ユーザー確認済み：一度3に引き上げたが、
# ちょうど2点共有している組も箱で強調したいとのことで2に戻した -- 2点
# 以上共有している組は全て箱で囲む)。INSET_PAD_FRACで外接矩形の外側に
# 少し余白を持たせる(点のマーカーぎりぎりに箱の辺が重ならないように)。
# badge_selection問わず(前版はbadge_selection="zoom_box"限定だったが、
# この分析自体はbadge選択方式と無関係なので2026-10-17に限定を外した)
# 常に実行する。
#
# ROWSには実際にannotateされた点のみを渡す(d_drawn/dd--描画されている
# 点全部--を試したところ、annotateされていない点まで箱に入るのはノイズ
# が多すぎるとユーザーからフィードバックがあり、2026-10-17に絞った)が、
# aggregate_locators=TRUEの時のbadge_rowsをそのまま渡してはいけない --
# .aggregate_locator_groups()が既に同じ(n_target, n_not_target, kind)
# キーで1行にまとめた後なので、ここでもう一度同じキーでグループ化しても
# 各グループが必ず1行になり、箱が一つも描けない(2026-10-17ユーザー報告
# で発覚："框完全消失了")。呼び出し元(draw_point_vs_target)は
# aggregate_locators=FALSEでselect_badgesをもう一度呼んだ
# tie_box_rows(「annotateされた点」ではあるが集約されていない、各点
# 自身の実座標を保ったままの行)を渡すこと -- badge_rows自体(集約後、
# badge表示用)は変更しない。
#
# .compute_x_tie_boxes()はrect()を直接描かず、{key: {x1,y1,x2,y2}}
# (dataの座標系のまま)を返すだけ -- .draw_badges()自身の連結線も、この
# 点のbadgeがどこかの箱に属している時は箱の辺から出発させたい
# (2026-10-17、ユーザー確認済み：「線がもし框なら、框のedgeから出発」)
# ので、計算結果を.draw_x_tie_boxes()(実際に描く)と.draw_badges()
# (連結線の出発点を決める)の両方で共有できるよう、計算と描画を分離した。
# POINT_RADIUS_IN(ROWSと同じ長さ、各点を実際に描いた時のpch記号の半径を
# インチで)はオプション -- 渡されなければ0扱い(旧来通り、usrのごく小さい
# 一定割合だけの余白)。渡す理由(2026-10-17ユーザー報告"放大模式下没有
# 和node一样被扩大，导致node溢出框体")：INSET_PAD_FRACはusrの範囲に対する
# 固定割合なので、インチ換算すると「パネル幅の一定%」という、ズーム
# レベルに関係なく常に同じ絶対サイズになる(usrが狭くなるほど1データ単位
# あたりのインチ数が増えるので、data単位の余白をインチに戻すと両者が
# 打ち消し合う)。一方render_zoomed_scatter_panel側の点のcexはzoom_factor
# (最大4倍)で別途拡大されるため、ズームパネルでは記号の実際の半径が
# パネル幅に対して不釣り合いに大きくなり、この固定インチ幅の余白では
# 包みきれなくなる。各グループの最大記号半径(呼び出し元がそのパネルで
# 実際に使ったcex式と同じ式で算出し、インチのまま渡す)をデータ単位に
# 換算して、INSET_PAD_FRAC由来の余白とのmaxを取ることで、記号がどれだけ
# 拡大されていても箱が必ず記号を覆うようにする。
.compute_x_tie_boxes <- function(rows, point_radius_in = NULL, inset_pad_frac = 0.01,
                                  min_group_size = 2) {
  if (nrow(rows) < min_group_size) return(list())
  usr <- par("usr"); pin <- par("pin")
  xpad <- diff(usr[1:2]) * inset_pad_frac
  ypad <- diff(usr[3:4]) * inset_pad_frac
  x_per_in <- diff(usr[1:2]) / pin[1]
  y_per_in <- diff(usr[3:4]) / pin[2]
  if (is.null(point_radius_in)) point_radius_in <- rep(0, nrow(rows))
  key <- paste(rows$n_target, rows$n_not_target, sep = "/")
  groups <- split(seq_len(nrow(rows)), key)
  boxes <- list()
  for (k in names(groups)) {
    idx <- groups[[k]]
    if (length(idx) < min_group_size) next
    x_rng <- range(rows$frequency[idx])
    y_rng <- range(rows$keyness[idx])
    if (diff(x_rng) <= 0 && diff(y_rng) <= 0) next
    max_radius_in <- max(point_radius_in[idx], 0)
    pad_x <- max(xpad, max_radius_in * x_per_in)
    pad_y <- max(ypad, max_radius_in * y_per_in)
    boxes[[k]] <- list(x1 = x_rng[1] - pad_x, y1 = y_rng[1] - pad_y,
                        x2 = x_rng[2] + pad_x, y2 = y_rng[2] + pad_y)
  }
  boxes
}

# cex(文字/記号サイズの倍率)から、その記号(pch 0:25系)がおおよそ占める
# 半径をインチで概算する -- Rのpch記号は1文字分の高さ(par("cin")[2])を
# 基準にcex倍されるので、直径をcex*cin[2]、半径をその半分とする近似。
# ピクセル単位で厳密に一致させる必要はない(.compute_x_tie_boxes側は
# これにinset_pad_fracの余白もさらに足すので、多少大きめに見積もる分には
# 問題なく、小さすぎる方だけが実害になる)。
.cex_to_radius_in <- function(cex) cex * par("cin")[2] * 0.5

# 単純な箱のみ描く(badge=文字ラベル付きではない)、グループを選別する
# 閾値/上限は.compute_x_tie_boxes側のMIN_GROUP_SIZE以外に設けない
# (ユーザー確認済み)。灰色の破線(2026-10-17、赤の実線から変更、
# ユーザー確認済み)。
.draw_x_tie_boxes <- function(boxes, border_col = "#999999") {
  for (box in boxes) {
    rect(box$x1, box$y1, box$x2, box$y2, border = border_col, lty = 2, lwd = 1, xpd = NA)
  }
  invisible(NULL)
}

# 源框(主图里，标出被放大的区域)、目标框(放大图本身的绘图区)四个角一一
# 连线——两个矩形分别活在两套不同的par("usr")坐标系里(主图 vs 放大图)，
# 没法直接用data坐标连线，要先换算成NDC(跟fig=/mar=/画布大小都无关、只跟
# 设备大小挂钩的统一坐标系，grconvertX/Y(..., to="ndc")就是干这个的)。
# 调用方要在各自还活跃(par("usr")对应哪张图)的时候把两组NDC坐标都算好、
# 存下来，等两张图都画完了才能调这个函数——最后叠一层
# fig=c(0,1,0,1)+mar=c(0,0,0,0)的空白绘图区，让plot.window(0:1,0:1)的
# user坐标正好等于NDC，才能在这层上直接画跨越两张图的连线。跟已有的
# neutral(lty=2)/median(lty=1)参考线用不同的虚线样式(lty="1343"，点划
# 混合)，避免看起来像同一套参考线。
.draw_projection_lines <- function(src_ndc_x, src_ndc_y, dest_ndc_x, dest_ndc_y, col) {
  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
  plot.new()
  # xaxs/yaxs必须是"i"(精确匹配，不自动扩展)——plot.window()默认的"r"
  # 会在给定范围两端各自动加4%留白，导致这层0~1的user坐标系其实比真正的
  # NDC(0,1)更宽，src/dest_ndc_*(已经是精确的NDC值)拿去在这里画点/线时
  # 会被这层多余的留白拉偏，跟两张图各自画的框对不上("目前的虚线和主图
  # 中框定的4个角没有接上，位置有些分离"，2026-09-15定位到这个根因)。
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  corners <- expand.grid(xi = 1:2, yi = 1:2)
  for (k in seq_len(nrow(corners))) {
    xi <- corners$xi[k]; yi <- corners$yi[k]
    segments(src_ndc_x[xi], src_ndc_y[yi], dest_ndc_x[xi], dest_ndc_y[yi],
             col = col, lty = "1343", lwd = 0.8, xpd = NA)
  }
}

# ---- inset下半：主图(不画badge，只标出被放大的区域) ----
render_inset_source_panel <- function(d, ranges, nofreq, axis_label, point_label, breaks,
                                       zbox, inset_border_col, value_label, not_value_label,
                                       colors, default_label, kld_variant = NULL) {
  .draw_axes(d, ranges, nofreq, axis_label, point_label, breaks$x_break, breaks$y_break,
             kld_variant = kld_variant)
  points(d$frequency, d$keyness, pch = d$point_pch,
         cex = .size_scale(d$node_size), col = d$point_col)
  .draw_source_box(zbox, inset_border_col)
  # 必须趁现在(还没切换par(fig=)去画上半之前)算好NDC，晚了的话下半的
  # par("usr")就不再是"当前"绘图区的坐标系了
  src_ndc_x <- grconvertX(c(zbox$x1, zbox$x2), from = "user", to = "ndc")
  src_ndc_y <- grconvertY(c(zbox$y1, zbox$y2), from = "user", to = "ndc")
  pc <- panel_coord_converters()
  .draw_axis_direction_arrows(value_label, not_value_label)
  .draw_size_legend(d$n_target, d$node_size, pc$pin, pc$in_to_data_x, pc$in_to_data_y,
                     pch = .representative_pch(d$point_pch))
  .draw_shape_legend(d, pc$pin, pc$in_to_data_x, pc$in_to_data_y)
  .draw_legend(colors$type_ramps, colors$signed_range, colors$legend_cols, default_label, pc$pin,
               pc$in_to_data_x, pc$in_to_data_y, colors$default_palette)
  list(src_ndc_x = src_ndc_x, src_ndc_y = src_ndc_y)
}

# ---- inset上半(single_locator时)：整块面板换成一张网络图 ----
# Use the same R renderer in native R and WebAssembly; no subprocesses.
render_network_glyph_panel <- function(badge_rows, inset_border_col, value_label, not_value_label,
                                        network_scale, network_node_gap) {
  members <- if (!is.null(badge_rows$group_size)) badge_rows$members[[1]] else badge_rows
  if (nrow(members) > 500) {
    strength <- pmax(abs(members$assoc_w1_to_w2), abs(members$assoc_w2_to_w1))
    members <- members[strength >= quantile(strength, 0.99, na.rm = TRUE), , drop = FALSE]
  }
  field1 <- if (is.null(members$field1)) NA_character_ else members$field1
  field2 <- if (is.null(members$field2)) NA_character_ else members$field2
  w1 <- .concept_display_text(members$w1, field1)
  w2 <- .concept_display_text(members$w2, field2)
  nodes <- unique(c(w1, w2))
  v12 <- members$assoc_w1_to_w2
  v21 <- members$assoc_w2_to_w1
  c1 <- ifelse(v12 > v21, "#111111", "#CC2200")
  c2 <- ifelse(v12 > v21, "#CC2200", "#111111")
  edges <- data.frame(from = c(w1, w2), to = c(w2, w1))
  graph <- igraph::graph_from_data_frame(edges, directed = TRUE,
                                        vertices = data.frame(name = nodes))
  # A circular seed makes the force-directed layout reproducible without
  # changing the caller's random-number generator state.
  seed <- igraph::layout_in_circle(graph)
  layout <- if (length(nodes) > 1) igraph::layout_with_fr(graph, coords = seed, niter = 200) else seed
  node_colors <- ifelse(nodes %in% c(w1[c1 == "#CC2200"], w2[c2 == "#CC2200"]),
                        "#CC2200", "#111111")
  fields <- c(rep_len(field1, length(w1)), rep_len(field2, length(w2)))[match(nodes, c(w1, w2))]
  shapes <- ifelse(is.na(fields) | fields == "entry", "circle", "square")
  strength <- abs(c(v12, v21))
  max_strength <- max(strength, na.rm = TRUE)
  widths <- if (is.finite(max_strength) && max_strength > 0) pmax(0.8, 4 * strength/max_strength) else rep(0.8, length(strength))
  igraph::plot.igraph(graph, layout = layout,
      vertex.label = nodes, vertex.label.family = getOption("kld.plot.font_family", "sans"),
      vertex.label.cex = 0.7, vertex.label.color = node_colors,
      vertex.color = "white", vertex.frame.color = node_colors, vertex.shape = shapes,
      vertex.size = (8 + 4 * network_scale) / max(0.2, network_node_gap),
      edge.color = c(c1, c2), edge.width = widths, edge.arrow.size = 0.25,
      edge.curved = 0.15, margin = 0.2)
  box(col = inset_border_col, lwd = 1.5)
  legend("topleft", legend = sprintf("%s:%d / %s:%d", value_label, badge_rows$n_target[1],
         not_value_label, badge_rows$n_not_target[1]), bty = "n", cex = 0.65)
  pc <- panel_coord_converters()
  list(pending = list(),
       dest_ndc_x = grconvertX(pc$usr[1:2], from = "user", to = "ndc"),
       dest_ndc_y = grconvertY(pc$usr[3:4], from = "user", to = "ndc"),
       in_to_data_x = pc$in_to_data_x, in_to_data_y = pc$in_to_data_y, pin = pc$pin)
}

# ---- inset上半(正常情况)：独立完整坐标轴+badge的放大图 ----
render_zoomed_scatter_panel <- function(d, badge_rows, tie_box_rows, zbox, ranges, nofreq, axis_label, point_label,
                                         inset_border_col, badge_selection, badge_layout,
                                         value_label, not_value_label, neutral, max_group_members,
                                         zoom_annotated_only = FALSE, repel_annotated_only = FALSE,
                                         kld_variant = NULL) {
  box_ranges <- list(xr = c(zbox$x1, zbox$x2), yr = c(zbox$y1, zbox$y2),
                      xpad = 0, ypad = 0, neutral = ranges$neutral)
  .draw_axes(d, box_ranges, nofreq, axis_label, point_label, NULL, NULL, axis_col = inset_border_col,
             kld_variant = kld_variant)
  dd <- d[d$frequency >= zbox$x1 & d$frequency <= zbox$x2 &
            d$keyness >= zbox$y1 & d$keyness <= zbox$y2, ]
  if (should_hide_unannotated_in_zoom(zoom_annotated_only) ||
        should_hide_unannotated_for_repel(repel_annotated_only, badge_layout)) {
    dd <- filter_to_annotated_points(dd, badge_rows)
  }
  # 点的cex要按x/y两个方向各自的放大倍数(主图实际显示跨度/zbox跨度)算一个
  # 综合缩放系数(sqrt(x*y)几何平均，跟.size_scale本身面积换算同一个习惯)，
  # 不然点的"直径"不会随框被放大而变大(用户反馈"扩大图中的点没有等比例扩大")。
  # 但badge候选点挤得越紧(取景框越窄)这个倍数越夸张，可以轻松到20~30倍，
  # 点会炸到不成比例的大小(2026-09-26用户反馈"否则会过于大")——用
  # .capped_zoom_factor()封顶，封顶后仍是整个面板共用的同一个数，跟
  # .size_scale()各点自己的大小比例无关，不会出现"有的点封顶有的没封顶"
  # 这种各点不一致的情况(2026-09-26讨论确认这个理解)。
  zoom_x <- (diff(ranges$xr) + 2 * ranges$xpad) / diff(c(zbox$x1, zbox$x2))
  zoom_y <- (diff(ranges$yr) + 2 * ranges$ypad) / diff(c(zbox$y1, zbox$y2))
  zoom_factor <- .capped_zoom_factor(sqrt(zoom_x * zoom_y))
  points(dd$frequency, dd$keyness, pch = dd$point_pch,
         cex = zoom_factor * .size_scale(dd$node_size), col = dd$point_col)
  box(col = inset_border_col, lwd = 1.5)

  pc <- panel_coord_converters()
  badge_cex <- 0.55
  fam <- getOption("kld.plot.font_family", "Hiragino Mincho ProN")
  label_gap_in <- strheight("あ", cex = badge_cex, family = fam) * 0.5
  # 投影线目标端连到放大图这边的框角——这边box就是整个坐标轴窗口，
  # par("usr")本身正好等于box，直接用pc$usr换算NDC即可
  dest_ndc_x <- grconvertX(pc$usr[1:2], from = "user", to = "ndc")
  dest_ndc_y <- grconvertY(pc$usr[3:4], from = "user", to = "ndc")

  tie_boxes <- .compute_x_tie_boxes(
    tie_box_rows,
    point_radius_in = .cex_to_radius_in(zoom_factor * .size_scale(tie_box_rows$node_size))
  )
  .draw_x_tie_boxes(tie_boxes)
  badge_colors <- .badge_color_from_point(badge_rows)
  .draw_badges(badge_rows, fill = badge_colors$fill, border = badge_colors$border_col,
               text_col = badge_colors$text_col, word_col = badge_colors$word_col,
               badge_layout = badge_layout, d_all = dd,
               value_label = value_label, not_value_label = not_value_label,
               badge_cex = badge_cex, fam = fam, label_gap_in = label_gap_in,
               x_per_in = pc$x_per_in, y_per_in = pc$y_per_in,
               data_to_in_x = pc$data_to_in_x, data_to_in_y = pc$data_to_in_y,
               in_to_data_x = pc$in_to_data_x, in_to_data_y = pc$in_to_data_y,
               pin = pc$pin, neutral = neutral, max_group_members = max_group_members,
               tie_boxes = tie_boxes)
  .draw_axis_direction_arrows(value_label, not_value_label)
  list(pending = list(), dest_ndc_x = dest_ndc_x, dest_ndc_y = dest_ndc_y,
       in_to_data_x = pc$in_to_data_x, in_to_data_y = pc$in_to_data_y, pin = pc$pin)
}

# 主图/放大图两块正方形panel各自的fig区间+画布该往哪个方向翻倍——
# width_multiplier=2时画布整体加宽(左右分栏)，height_multiplier=2时画布
# 整体加高(上下叠放)，两者互斥(每次只翻一个方向，另一个方向保持1倍，
# 两块panel永远是正方形，不会因为翻倍方向不同就被拉伸变形)。"top"是
# 2026-09-21就有的默认布局("扩大图改为放在上方")，width_multiplier这个
# 参数当时就是为了给"left/right"这类以后可能要支持的方向预留的(见
# .open_canvas自己的注释"width_multiplier保留但目前没有调用方在用")，
# 2026-09-29正式用上，加了"bottom"/"left"/"right"三个新方向(用户要求
# "可以调整放大框放在主图的上下左右")。
.inset_fig_regions <- function(inset_position) {
  switch(inset_position,
    top = list(main = c(0, 1, 0, 0.5), inset = c(0, 1, 0.5, 1),
               width_multiplier = 1, height_multiplier = 2),
    bottom = list(main = c(0, 1, 0.5, 1), inset = c(0, 1, 0, 0.5),
                  width_multiplier = 1, height_multiplier = 2),
    left = list(main = c(0.5, 1, 0, 1), inset = c(0, 0.5, 0, 1),
                width_multiplier = 2, height_multiplier = 1),
    right = list(main = c(0, 0.5, 0, 1), inset = c(0.5, 1, 0, 1),
                 width_multiplier = 2, height_multiplier = 1),
    stop("inset_position must be one of top/bottom/left/right, got: ", inset_position)
  )
}

# inset=TRUE的完整渲染路径，跟非inset路径(draw_point_vs_target主体)完全
# 独立，互不影响。.inset_fig_regions()按inset_position算出画布该往哪个
# 方向翻倍、主图/放大图各自的fig区间，分别调用.draw_axes——每一半单独看，
# 跟一张普通的(非inset)图用同一套margin/字号，视觉上完全一致。校准标记/
# 投影线/网络图svg拼接这几处都只读"当前panel自己的坐标系"(par("usr")/
# par("pin")这类)，不假设主图在下还是放大图在上，inset_position换方向
# 不需要改这几处。
.render_with_inset <- function(d, badge_rows, tie_box_rows, colors, output, ranges, breaks, neutral, nofreq,
                                 axis_label, point_label, value_label, not_value_label,
                                 default_label, badge_selection, badge_layout,
                                 canvas_size, inset_border_col, inset_pad_frac,
                                 inset_position = "top",
                                 zoom_box = NULL, max_group_members = 50, network_scale = 0.01,
                                 network_node_gap = 1,
                                 zoom_annotated_only = FALSE, repel_annotated_only = FALSE,
                                 kld_variant = NULL) {
  # colors/d$point_col由调用方(draw_point_vs_target)在选badge_rows之前
  # 算好、传进来——不在这里重算，badge_rows本身就是d的一份快照，如果这里
  # 重算一次d$point_col，也补救不了badge_rows早就没有这一列的事实
  # (2026-09-15修的bug，见调用方那边的注释)。
  # 局部变量不能叫box——base R自己有个同名的box()函数(画四条边的框，
  # 下面上半就用得到)，叫box会把那个函数遮蔽掉
  # ZOOM_BOX给了就直接用调用方指定的取景范围(比如"y>1,x>0"这种手动划定
  # 的区域，2026-09-19需求)，不再从badge_rows的外接框自动算——手动指定的
  # 范围本身已经是完整区间，不需要.compute_zoom_box那套"退化维度找最小
  # 宽度兜底"的处理(那是自动算badge外接框时才会遇到的问题)。
  # badge_rows(aggregate_locators=TRUEの時は集約済み)ではなく、必ず
  # 集約前のtie_box_rowsを渡す -- 同じ(n_target, n_not_target, kind)に
  # 束ねられたグループは代表行1つの座標しか持たないため、badge_rowsを
  # 渡すとグループ内の他メンバー(色は塗られているが代表行より低い
  # dispersionを持つ点)がこの外接矩形の計算から漏れる(2026-10-19
  # ユーザー報告:"有些低dispersion的词被标红了，当时没有被框起来")。
  # tie_box_rowsは.compute_x_tie_boxes向けに同じ理由で既に非集約のまま
  # 渡されている(上のコメント参照) -- ここも同じ原則に揃える。
  zbox <- if (!is.null(zoom_box)) zoom_box else .compute_zoom_box(tie_box_rows, ranges, inset_pad_frac)

  # 扩大图范围里所有候选badge点(不管聚不聚合、不管是不是同一组)汇总起来
  # 是不是只对应一个locator(n_target/n_not_target完全相同，且是同一种pair
  # 组合field1+field2——不同种类不算同一个locator，见.aggregate_locator_
  # groups()同一处注释，2026-09-23"badge 不同种类的pair和condition保持不
  # 聚合")——这种极端情况下"多个badge摆在坐标轴上各自找空位"已经没有意义
  # (反正只有一个点、一个坐标)，改成整个上半面板直接画一张网络图，不画
  # 坐标轴、只留一个框(用户2026-09-22要求)。aggregate_locators=TRUE时同一
  # locator+同一种类的多行早就被.aggregate_locator_groups()合并成1行了，
  # 这时nrow(badge_rows)==1本身就等价于这个判断；FALSE时badge_rows可能还是
  # 好几行(都是同一个locator+种类没合并)，所以统一按locator+kind去重之后
  # 数unique值，不直接判nrow==1。
  is_pair <- !is.null(badge_rows$w1)
  locator_kind_key <- function(rows) {
    kind <- if (!is.null(rows$field1)) paste(rows$field1, rows$field2, sep = "|") else ""
    paste(rows$n_target, rows$n_not_target, kind, sep = "\t")
  }
  single_locator <- is_pair && nrow(badge_rows) > 0 &&
    length(unique(locator_kind_key(badge_rows))) == 1

  regions <- .inset_fig_regions(inset_position)
  .open_canvas(output, canvas_size, width_multiplier = regions$width_multiplier,
               height_multiplier = regions$height_multiplier)
  # 跟draw_point_vs_target()非inset路径同一个理由：正常情况由下面显式的
  # dev.off()关闭设备，splice_ok只是避免on.exit兜底再关一次已经关掉的设备
  splice_ok <- FALSE
  on.exit(if (!splice_ok) dev.off())

  # 两块正方形panel各占设备翻倍那个方向总长度的一半，具体谁在哪一半由
  # regions$main/regions$inset(即inset_position)决定。
  # 主图：不画badge，只画虚线框标出被放大的区域
  par(fig = regions$main)
  source_panel <- render_inset_source_panel(d, ranges, nofreq, axis_label, point_label, breaks,
                                             zbox, inset_border_col, value_label, not_value_label,
                                             colors, default_label, kld_variant = kld_variant)

  # 放大图：正常情况下是独立的完整坐标轴+badge；single_locator时整块
  # 面板换成一张网络图(用户2026-09-22要求"这个网络图不要坐标轴，只需要一个
  # 框")。坐标轴/框都用inset_border_col(红)而不是默认黑色(2026-09-15确认
  # "轴上的标度改为红色")。
  par(fig = regions$inset, new = TRUE)

  panel <- if (single_locator) {
    render_network_glyph_panel(badge_rows, inset_border_col, value_label, not_value_label,
                                network_scale, network_node_gap)
  } else {
    render_zoomed_scatter_panel(d, badge_rows, tie_box_rows, zbox, ranges, nofreq, axis_label, point_label,
                                 inset_border_col, badge_selection, badge_layout,
                                 value_label, not_value_label, neutral, max_group_members,
                                 zoom_annotated_only = zoom_annotated_only,
                                 repel_annotated_only = repel_annotated_only,
                                 kld_variant = kld_variant)
  }
  pending <- panel$pending

  # 校准标记要在.draw_projection_lines()之前画——那个函数会把par(fig=)切到
  # 覆盖整个设备的[0,1]全画布叠加层(画src/dest两个框四角的连线用)，画完
  # 之后"当前绘图区"就不再是上半这块面板了，panel$in_to_data_x/y这两个闭包
  # 虽然算出来的还是上半坐标系下的正确数值，但rect()会拿这个数值去套用
  # "当前"那个[0,1]叠加层的坐标系，位置整个错掉(2026-09-21用真实数据跑出来
  # 的bug："找不到svg校准标记")。校准标记必须在上半坐标系还是"当前"绘图区
  # 时画完。
  marks_in <- if (length(pending) > 0) {
    .draw_svg_calibration_marks(panel$in_to_data_x, panel$in_to_data_y, panel$pin)
  } else NULL

  .draw_projection_lines(source_panel$src_ndc_x, source_panel$src_ndc_y,
                          panel$dest_ndc_x, panel$dest_ndc_y, inset_border_col)

  dev.off()
  splice_ok <- TRUE
  if (length(pending) > 0) {
    .splice_network_glyphs(output, pending, marks_in)
  }
}

.kv_pairs <- function(tokens) {
  kv <- strsplit(tokens, "=", fixed = TRUE)
  setNames(vapply(kv, `[`, character(1), 2), vapply(kv, `[`, character(1), 1))
}

.args_from_cli <- function() {
  .kv_pairs(commandArgs(trailingOnly = TRUE))
}

# Pure: ARGS(一份key=value的具名字符向量，跟.args_from_cli()的输出同一个
# 形状) -> draw_point_vs_target()自己需要的完整kwargs list，带跟CLI入口
# (下面sys.nframe()==0那段)完全一致的类型转换/默认值。从那段代码里拆出来
# 单独成一个函数，是为了让app.R的"直接粘贴cmdline"输入模式(2026-10-17
# 新增)能调用同一份映射逻辑，而不是另外抄一份容易跟真正的CLI契约
# 走着走着就不一致——跟这个文件里ALL_PARAM_KEYS"两处共用一份列表"同一个
# 道理(app.R自己注释里也有这段解释)。OUTPUT单独传是因为两个调用方对它
# 的处理不一样：CLI入口直接读args[["output"]]；app.R的cmdline模式会
# 用自己的tempfile覆盖掉(不管用户粘贴的命令里写了什么output=)。
.draw_point_vs_target_args <- function(args, output = args[["output"]]) {
  .arg <- function(name, default = NULL) {
    if (name %in% names(args)) args[[name]] else default
  }
  list(
    input = args[["input"]], output = output,
    point_label = .arg("point_label", "point"),
    axis_label = .arg("axis_label", "target"),
    value_label = .arg("value_label", "value"),
    not_value_label = .arg("not_value_label", paste0("¬", .arg("value_label", "value"))),
    nofreq = tolower(.arg("nofreq", "false")) == "true",
    badge_layout = .arg("badge_layout", "empty_region"),
    annotate_words_file = .arg("annotate_words_file"),
    highlight_file = .arg("highlight_file"),
    highlight_pairs_file = .arg("highlight_pairs_file"),
    annotate_pairs_file = .arg("annotate_pairs_file"),
    annotate_nontarget_only = tolower(.arg("annotate_nontarget_only", "false")) == "true",
    annotate_point_exclude = if (!is.null(.arg("annotate_point_exclude"))) {
      strsplit(.arg("annotate_point_exclude"), ",", fixed = TRUE)[[1]]
    } else NULL,
    highlight_frontier = tolower(.arg("highlight_frontier", "false")) == "true",
    default_label = .arg("default_label", .arg("point_label", "point")),
    badge_selection = .arg("badge_selection", "words"),
    frontier_n = as.numeric(.arg("frontier_n", "15")),
    inset = tolower(.arg("inset", "false")) == "true",
    inset_position = .arg("inset_position", "top"),
    inset_pad_frac = as.numeric(.arg("inset_pad_frac", "0.15")),
    zoom_annotated_only = tolower(.arg("zoom_annotated_only", "false")) == "true",
    repel_annotated_only = tolower(.arg("repel_annotated_only", "false")) == "true",
    zoom_x1 = .arg("zoom_x1"),
    zoom_x2 = .arg("zoom_x2"),
    zoom_y1 = .arg("zoom_y1"),
    zoom_y2 = .arg("zoom_y2"),
    edge_threshold = .arg("edge_threshold"),
    dispersion_threshold = .arg("dispersion_threshold"),
    node_dispersion_threshold = .arg("node_dispersion_threshold"),
    frequency_threshold = .arg("frequency_threshold"),
    node_frequency_threshold = .arg("node_frequency_threshold"),
    aggregate_locators = tolower(.arg("aggregate_locators", "true")) == "true",
    max_group_members = as.numeric(.arg("max_group_members", "50")),
    network_scale = as.numeric(.arg("network_scale", "0.01")),
    network_node_gap = as.numeric(.arg("network_node_gap", "1")),
    point_field = .arg("point_field", "entry"),
    kld_variant = .arg("kld_variant"),
    frequency_type = .arg("frequency_type")
  )
}

# 只有直接被Rscript调用(不是被source())时才跑CLI入口 -- run_meta_config_plots.R
# source()这个文件是为了复用draw_point_vs_target()，不希望连带触发这一段。
if (sys.nframe() == 0 && !interactive()) {
  do.call(draw_point_vs_target, .draw_point_vs_target_args(.args_from_cli()))
}
