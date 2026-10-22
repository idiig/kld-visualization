FROM rocker/r-ver:4.5.2 AS runtime
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \
    libfreetype6-dev libpng-dev libtiff-dev libjpeg-dev libharfbuzz-dev \
    libfribidi-dev libglpk-dev libarchive-dev fonts-noto-cjk \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY install.R ./
RUN Rscript install.R
COPY . .
ENV KLD_PLOT_FONT="Noto Serif CJK JP"
EXPOSE 3838
CMD ["Rscript", "-e", "shiny::runApp('.', host='0.0.0.0', port=3838, launch.browser=FALSE)"]

FROM runtime AS pages-build
RUN Rscript scripts/build_site.R

FROM nginx:alpine AS pages
COPY --from=pages-build /app/_site /usr/share/nginx/html
EXPOSE 80
