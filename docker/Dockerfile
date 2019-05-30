From debian:stretch

RUN apt-get update && \
  apt-get install -y build-essential \
    wget \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libncurses-dev \
    libbz2-dev \
    liblzma-dev \
    libxml2-dev \
    libgfortran3 \
    libcairo2-dev \
    zip \
    python3-dev \
    python3-pip \
    pandoc \
    r-base \
    git

# Install all software under /opt/software:
RUN mkdir -p /opt/software 

# Get the miniconda installer script and run it in silent mode:
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt/software/miniconda.sh && \
  bash /opt/software/miniconda.sh -b -p /opt/software/miniconda
ENV PATH="/opt/software/miniconda/bin:${PATH}"

# install kallisto
RUN conda config --add channels defaults && \
  conda config --add channels bioconda && \
  conda config --add channels conda-forge && \
  conda install -y kallisto && \
  conda install -y r-sleuth

# Install fastQC:
RUN cd /opt/software && \
  wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip && \
  unzip fastqc_v0.11.8.zip && \
  chmod +x /opt/software/FastQC/fastqc
ENV PATH="/opt/software/FastQC:${PATH}"

# Install libStatGen and fastQValidator:
RUN cd /opt/software && \
  git clone https://github.com/statgen/libStatGen.git && \
  cd libStatGen && \
  make && \
  cd .. && \
  git clone https://github.com/statgen/fastQValidator.git && \
  cd fastQValidator && \
  make && \
  make install

# Install DESeq2 and dependencies:
RUN R -e "install.packages('ggplot2', repos = 'http://cran.us.r-project.org');"
#RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('RCurl'); biocLite('XML');"
#RUN R -e "install.packages('reshape2', repos = 'http://cran.us.r-project.org'); install.packages('tidyverse', repos = 'http://cran.us.r-project.org'); install.packages('ggdendro', repos = 'http://cran.us.r-project.org'); install.packages('gplots', repos = 'http://cran.us.r-project.org'); install.packages('ggplot2', repos = 'http://cran.us.r-project.org');"
#ADD deseq2.R /opt/software/
#ADD make_dge_plots.py /opt/software/
#ADD make_figures.R /opt/software/
#ADD draw_heatmap.R /opt/software/
#ADD draw_pca.R /opt/software/
#ADD contrast_independent_figures.R /opt/software/

# Install some Python3 libraries:
ADD requirements.txt /opt/software/
RUN pip3 install -r /opt/software/requirements.txt

# Install multiqc:
RUN pip3 install multiqc

# Install Jinja2 for templating.  Likely installed with other python packages, but just to be sure:
RUN pip3 install Jinja2

# for multiqc to read commandline:
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

# The scripts for performing pre-workflow checks:
RUN mkdir /opt/software/precheck
ADD perform_precheck.py /opt/software/precheck/
ADD check_fastq.py /opt/software/precheck/

# The script for generating the markdown report:
ADD generate_report.py /usr/local/bin/
RUN chmod +x /usr/local/bin/generate_report.py

# Add the report template and styling:
RUN mkdir /opt/report
ADD report.md /opt/report/
ADD report.css /opt/report/

ENTRYPOINT ["/bin/bash"]

