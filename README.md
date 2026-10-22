# KLD visualization

在线绘图：**https://idiig.github.io/kld-visualization/**

源码：https://github.com/idiig/kld-visualization

独立的 R 绘图和在线调参 App。GitHub Pages 使用 Shinylive/webR 在浏览器运行 R；本地与 Docker 使用普通 Shiny。两种模式共用绘图代码。

仓库只包含展示代码、已有绘图结果数据和注释列表，不包含统计计算管线、原始语料、索引或计算配置。

## Nix

```sh
nix develop
Rscript tests/smoke.R
Rscript -e 'shiny::runApp(".", host="127.0.0.1", port=3838)'
```

打开 http://localhost:3838 。使用不启用 flakes 的 Nix 时，在 `nix` 后加 `--extra-experimental-features 'nix-command flakes'`。

构建、预览 GitHub Pages 相同的浏览器版本：

```sh
nix develop --command Rscript scripts/build_site.R
nix develop --command Rscript -e 'httpuv::runStaticServer("_site", port=8080)'
```

打开 http://localhost:8080 。第一次构建需要联网下载 Shinylive 和 WebAssembly 包。

## Docker

```sh
docker compose up --build app
```

打开 http://localhost:3838 。直接构建和预览浏览器版：

```sh
docker compose --profile pages up --build pages
```

打开 http://localhost:8080 。Docker 和 Nix 是本地/构建环境；GitHub Pages 线上不运行容器或 R 服务。

## GitHub Pages

1. 将此目录作为独立 Git 仓库推送到 GitHub 的 `main` 分支。
2. 在仓库 **Settings → Pages → Build and deployment → Source** 选择 **GitHub Actions**。
3. `.github/workflows/pages.yaml` 自动执行绘图检查、导出 `_site`，然后发布。PR 只构建检查。

运行后的网址由 Actions 部署输出提供。无需填写固定仓库路径，支持项目 Pages 的子路径。

## 数据和绘图

`data/processed/visualization/*.rds` 是预计算 TSV 的无损序列化版本，保留全部行和列；`manifest.json` 记录行数、列名及原 TSV 的 MD5。App 中保留原 TSV 名称，已有命令仍可读取对应 RDS：

```sh
Rscript scripts/visualization/plot_point_vs_target.R \
  input=data/processed/visualization/entry-vs-gender-kokin.tsv \
  output=artifacts/example.svg kld_variant=norm frequency_type=term
```

仅重新打包外部已算好的 TSV（不执行统计计算）：

```sh
Rscript scripts/package_data.R /path/to/precomputed/visualization
```

`data/reference/` 为高亮和注释列表。自动标签是 App 内独立的展示元数据。

浏览器版保留参数调整、筛选、放大图、网络图、结果表、复现命令和 SVG 下载。网络图统一由 R/igraph 绘制，布局与原 Graphviz 版本不同；无需外部 `neato`。日文 PNG/SVG 使用仓库附带的 Noto Serif CJK JP 字体，许可证见 `fonts/OFL.txt`。本地环境的字体可通过 `KLD_PLOT_FONT` 配置。

浏览器按需下载所选表，已下载的表在当前浏览器会话中保留。大型词对表会占用较多内存，建议桌面浏览器。绘图由「描画」按钮触发。

## 浏览器验证

原生绘图检查为 `Rscript tests/smoke.R`。端到端检查使用 Chromium：

```sh
npm install --no-save playwright
npx playwright install chromium
KLD_SERVE=1 node tests/browser.mjs http://127.0.0.1:8080/
```

`KLD_SERVE=1` 自动启动/关闭 `_site` 的本地静态服务器，需要 Python 3；也可以先自行托管 `_site`，再向测试脚本传入网址。测试覆盖启动、数据选择、PNG 绘图、参数重绘、SVG 下载、结果表及词对绘图。
