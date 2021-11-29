FROM rocker/shiny:4.1.1

RUN apt-get update -y && apt-get install -y libssl-dev python3-pip libmagick++-dev libharfbuzz-dev libfribidi-dev libmariadbclient-dev libwebp-dev cargo libhdf5-dev libnetcdf-dev libavfilter-dev
RUN pip3 install statsmodels missingpy sklearn==0.21.3 pandas RegImpute tkinter
RUN R -e "install.packages('BiocManager', dependencies=TRUE, repos='https://cran.us.r-project.org')"
RUN R -e "BiocManager::install(c('MSnbase','limma','ROTS'))"
RUN R -e "install.packages(c('reshape','lattice','ggplot2','limma','ROTS','stringr','seqinr','ggrepel','MKmisc','gplots','gtools','magick','dplyr','shinythemes','mwshiny'), dependencies=TRUE, repos='https://cran.us.r-project.org')"

WORKDIR /srv/shiny-server/PEA

RUN chown -R shiny:shiny .

CMD ["/init"]

