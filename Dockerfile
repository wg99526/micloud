FROM openanalytics/r-base
#FROM wg99526/micloudtest

#MAINTAINER Tobias Verbeke "tobias.verbeke@openanalytics.eu"
MAINTAINER Won Gu "rndnjs526@gmail.com"


# system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libssl1.1 \
    libxml2-dev \
    build-essential \
    r-base-dev \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# system library dependency for the euler app
RUN apt-get update && apt-get install -y \
    libmpfr-dev \
    && rm -rf /var/lib/apt/lists/*

# basic shiny functionality
RUN R -e "install.packages(c('shiny', 'rmarkdown'), repos='https://cloud.r-project.org/')"

# install dependencies of the euler app
#RUN R -e "install.packages('Rmpfr', repos='https://cloud.r-project.org/')"

# install dependencies of the MiCloud app
RUN R -e "install.packages(c('seqinr', 'shinydashboard', 'dashboardthemes', 'tidyverse', 'plotly', 'shinyWidgets', 'shinyjs', 'googleVis', 'xtable'), repos='https://cloud.r-project.org/')"
RUN R -e "install.packages(c('DT', 'htmltools', 'phangorn', 'bios2mds', 'zip', 'zCompositions', 'dplyr', 'forestplot', 'quantreg', 'fossil', 'picante' ), repos='https://cloud.r-project.org/')"
RUN R -e "install.packages(c('entropart', 'lme4', 'lmerTest', 'broom.mixed', 'gee', 'geepack', 'dirmult', 'robustbase', 'robCompositions', 'BiasedUrn'), repos='https://cloud.r-project.org/')"
RUN R -e "install.packages(c('CompQuadForm', 'GUniFrac', 'ecodist', 'MiRKAT', 'gridExtra', 'ggplot2', 'patchwork', 'ggthemes', 'erer', 'DiagrammeR', 'stringr'), repos='https://cloud.r-project.org/')"
RUN R -e "install.packages(c('devtools', 'betareg', 'reticulate', 'nlme', 'glmmTMB', 'glmm', 'remotes'), repos='https://cloud.r-project.org/')"


#RUN R -e "install.packages(c('phyloseq', 'biomformat'), repos='https://bioconductor.org/biocLite.R')"
RUN R -e "remotes::install_github('joey711/phyloseq')"
RUN R -e "remotes::install_github('joey711/biomformat')"
#RUN R -e "remotes::install_github('GLMMMiRKAT', repos='https://github.com/hk1785/GLMM-MiRKAT.git')"
RUN R -e "remotes::install_github('hk1785/GLMM-MiRKAT')"
RUN R -e "remotes::install_github('nyiuab/NBZIMM')"
RUN R -e "install.packages(c('gridGraphics', 'compositions'), repos = 'https://cloud.r-project.org/')"

# copy the app to the image
#RUN mkdir /root/euler
#COPY euler /root/euler
RUN mkdir /root/app
COPY app /root/app
COPY Rprofile.site /usr/lib/R/etc/

#COPY app/sub.biom.Rdata /root/app

COPY app/Data/BMI/biom.MZ.BMI.Rdata /root/app
COPY app/Data/BMI/otu.tab.txt /root/app
COPY app/Data/BMI/sam.dat.txt /root/app
COPY app/Data/BMI/tax.tab.txt /root/app
COPY app/Data/BMI/tree.tre /root/app

COPY app/Data/val_physeq.rds /root/app

COPY app/MiDataProc.Data.Upload.R /root/app
COPY app/MiDataProc.Alpha.Cross.Sectional.R /root/app
COPY app/MiDataProc.Alpha.Longitudinal.R /root/app
COPY app/MiDataProc.Beta.Cross.Sectional.R /root/app
COPY app/MiDataProc.Beta.Longitudinal.R /root/app
COPY app/MiDataProc.Taxa.Cross.Sectional.R /root/app
COPY app/MiDataProc.Taxa.Longitudinal.R /root/app

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/root/app')"]
