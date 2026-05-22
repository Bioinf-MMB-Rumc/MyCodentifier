FROM broadinstitute/gatk:4.1.4.1

# Install miniconda to /miniconda
RUN curl -kLO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b && \
    rm Miniconda3-latest-Linux-x86_64.sh

# set the environment
ENV PATH="/miniconda/bin:$PATH" LC_ALL=C

# Install miniconda environment
RUN conda install -c conda-forge mamba=1.1.0 && \
    mamba update --all && \
    mamba create -n mycodentifier python=3.6.8 pysam=0.15.2 nextflow=22.10.4 simplejson=3.17.0 pandas=1.0.1 matplotlib=3.1.2 -c bioconda -c conda-forge

## activate by bash mycoprofiler env
RUN echo "source activate mycodentifier" > ~/.bashrc
ENV PATH=/miniconda/envs/mycodentifier/bin:${PATH}

# set shell
SHELL ["/bin/bash", "-c"]

# add codebase to docker
ADD ./bin/snpit-master /workflow/bin/snpit-master

# adding snpit data from git and
RUN cd /workflow/bin/snpit-master && \
    python setup.py install && \
    python setup.py test

#install nano
RUN apt-get update && \
    apt-get install -y nano && \
    apt-get install -y fontconfig
    
#unset JAVA_HOME
ENV NXF_JAVA_HOME=/miniconda/envs/mycodentifier