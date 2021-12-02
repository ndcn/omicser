FROM rocker/r-ver:4.1.2
LABEL maintainer="Rico Derks" \
      email="r.j.e.derks@lumc.nl"

# install some depencies
RUN apt-get update && apt-get install -y  git-core \
    libcurl4-openssl-dev \
    libgit2-dev libicu-dev \
    libpng-dev \
    libssl-dev \
    libxml2-dev \
    make \
    pandoc \
    pandoc-citeproc \
    python zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

# set CRAN repo
RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" >> /usr/local/lib/R/etc/Rprofile.site

# install R packages
RUN R -e 'install.packages(c("remotes", "reticulate"))'
RUN Rscript -e 'remotes::install_version("magrittr",upgrade="never", version = "2.0.1")'
RUN Rscript -e 'remotes::install_version("shiny",upgrade="never", version = "1.7.1")'
RUN Rscript -e 'remotes::install_version("stringr",upgrade="never", version = "1.4.0")'
RUN Rscript -e 'remotes::install_version("knitr",upgrade="never", version = "1.36")'
RUN Rscript -e 'remotes::install_version("dplyr",upgrade="never", version = "1.0.7")'
RUN Rscript -e 'remotes::install_version("waiter",upgrade="never", version = "0.2.4")'
RUN Rscript -e 'remotes::install_version("data.table",upgrade="never", version = "1.14.2")'
RUN Rscript -e 'remotes::install_version("tidyr",upgrade="never", version = "1.1.4")'
RUN Rscript -e 'remotes::install_version("shinyjs",upgrade="never", version = "2.0.0")'
RUN Rscript -e 'remotes::install_version("shinydashboardPlus",upgrade="never", version = "2.0.3")'
RUN Rscript -e 'remotes::install_version("rmarkdown",upgrade="never", version = "2.11")'
RUN Rscript -e 'remotes::install_version("plotly",upgrade="never", version = "4.10.0")'
RUN Rscript -e 'remotes::install_version("markdown",upgrade="never", version = "1.1")'
RUN Rscript -e 'remotes::install_version("golem",upgrade="never", version = "0.3.1")'
RUN Rscript -e 'remotes::install_version("ggdendro",upgrade="never", version = "0.1.22")'
RUN Rscript -e 'remotes::install_version("DT",upgrade="never", version = "0.19")'
RUN Rscript -e 'remotes::install_version("configr",upgrade="never", version = "0.3.5")'
RUN Rscript -e 'remotes::install_version("broom",upgrade="never", version = "0.7.10")'
RUN Rscript -e 'remotes::install_version("anndata",upgrade="never", version = "0.7.5.3")'
RUN Rscript -e 'remotes::install_github("jokergoo/ComplexHeatmap@6c3ae544adfbf19566c70546c9235cfb9efe5481")'
RUN Rscript -e 'remotes::install_version("shinyFiles",upgrade="never", version = "0.9.1")'

### install omicser package
# make temporary work directory
RUN mkdir /build_zone
# add everything to the container
ADD . /build_zone
WORKDIR /build_zone
# install omicser
RUN R -e 'remotes::install_local(upgrade="never")'
# install miniconda stuff etc.
RUN R -e 'source("docker_stuff/install_py_env.R")'
# copy the .Renviron file
RUN cp ./docker_stuff/.Renviron /root/.Renviron
# copy a app config file
RUN cp ./docker_stuff/app_config.yml /root/app_config.yml
# cleanup after myself
RUN rm -rf /build_zone

# go back to the home directory
WORKDIR /root

# what port to expose
EXPOSE 3939

# fire the app on startup
CMD R -e "options('shiny.port'=3939,shiny.host='0.0.0.0');omicser::run_docker()"
