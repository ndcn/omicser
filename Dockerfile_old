FROM rocker/r-ver:4.1.0
RUN apt-get update && apt-get install -y  git-core imagemagick libcairo2-dev libcurl4-openssl-dev libgit2-dev libicu-dev libmagic-dev libmagick++-dev libpng-dev libssl-dev libxml2-dev make pandoc pandoc-citeproc python zlib1g-dev && \rm -rf /var/lib/apt/lists/*
RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" >> /usr/local/lib/R/etc/Rprofile.site
RUN R -e 'install.packages("remotes")'
RUN Rscript -e 'remotes::install_version("magrittr",upgrade="never", version = "2.0.1")'
RUN Rscript -e 'remotes::install_version("stringr",upgrade="never", version = "1.4.0")'
RUN Rscript -e 'remotes::install_version("knitr",upgrade="never", version = "1.36")'
RUN Rscript -e 'remotes::install_version("shiny",upgrade="never", version = "1.7.1")'
RUN Rscript -e 'remotes::install_version("dplyr",upgrade="never", version = "1.0.7")'
RUN Rscript -e 'remotes::install_version("tidyr",upgrade="never", version = "1.1.4")'
RUN Rscript -e 'remotes::install_version("ggplot2",upgrade="never", version = "3.3.5")'
RUN Rscript -e 'remotes::install_version("waiter",upgrade="never", version = "0.2.4")'
RUN Rscript -e 'remotes::install_version("shinydashboard",upgrade="never", version = "0.7.2")'
RUN Rscript -e 'remotes::install_version("data.table",upgrade="never", version = "1.14.2")'
RUN Rscript -e 'remotes::install_version("thinkr",upgrade="never", version = "0.15")'
RUN Rscript -e 'remotes::install_version("shinyWidgets",upgrade="never", version = "0.6.2")'
RUN Rscript -e 'remotes::install_version("shinyjs",upgrade="never", version = "2.0.0")'
RUN Rscript -e 'remotes::install_version("shinydashboardPlus",upgrade="never", version = "2.0.3")'
RUN Rscript -e 'remotes::install_version("rmarkdown",upgrade="never", version = "2.11")'
RUN Rscript -e 'remotes::install_version("plotly",upgrade="never", version = "4.10.0")'
RUN Rscript -e 'remotes::install_version("markdown",upgrade="never", version = "1.1")'
RUN Rscript -e 'remotes::install_version("magick",upgrade="never", version = "2.7.3")'
RUN Rscript -e 'remotes::install_version("gridExtra",upgrade="never", version = "2.3")'
RUN Rscript -e 'remotes::install_version("golem",upgrade="never", version = "0.3.1")'
RUN Rscript -e 'remotes::install_version("ggdendro",upgrade="never", version = "0.1.22")'
RUN Rscript -e 'remotes::install_version("DT",upgrade="never", version = "0.19")'
RUN Rscript -e 'remotes::install_version("configr",upgrade="never", version = "0.3.5")'
RUN Rscript -e 'remotes::install_version("broom",upgrade="never", version = "0.7.9")'
RUN Rscript -e 'remotes::install_version("BiocManager",upgrade="never", version = "1.30.16")'
RUN Rscript -e 'remotes::install_version("anndata",upgrade="never", version = "0.7.5.3")'
RUN Rscript -e 'remotes::install_bioc("3.13/SingleCellExperiment",upgrade="never")'
RUN Rscript -e 'remotes::install_bioc("3.13/LoomExperiment",upgrade="never")'
RUN Rscript -e 'remotes::install_github("cellgeni/sceasy@f8f0628a280e0880ea94b00100b463e1f6ba1994")'
#RUN Rscript -e 'remotes::install_github("jokergoo/ComplexHeatmap@9c277dda153b5f9bfdeb76eb1b087cf81b08bf80")'
RUN Rscript -e 'remotes::install_bioc("3.13/ComplexHeatmap",upgrade="never")'
RUN mkdir /build_zone
ADD . /build_zone
WORKDIR /build_zone
RUN R -e 'remotes::install_local()'
RUN rm -rf /build_zone
EXPOSE 80
CMD R -e "options('shiny.port'=80,shiny.host='0.0.0.0');omicser::run_app()"
