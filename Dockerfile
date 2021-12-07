FROM rocker/shiny-verse:4.1.1

RUN apt-get update -y && apt-get install -y --no-install-recommends \
  cython3 \
  imagemagick \
  libmagick++-dev \
  libnetcdf-dev \
  python3-dev \
  python3-pip \
  python3-statsmodels \
  && rm -rf /var/lib/apt/lists/*

# RegImpute requires scikit-learn < 0.22.0

RUN pip3 install --no-cache-dir \
  missingpy==0.2.0 \
  RegImpute==0.0.2 \
  scikit-learn==0.21.3

RUN Rscript -e "BiocManager::install(c('MSnbase','limma','ROTS'), Ncpus=parallel::detectCores())"

RUN install2.r --error --skipinstalled \
  ggrepel \
  gplots \
  gtools \
  magick \
  MKmisc \
  mwshiny \
  reshape \
  seqinr \
  shinythemes

WORKDIR /srv/shiny-server/PEA

COPY --chown=shiny:shiny . .
