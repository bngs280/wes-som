# Use an official Python runtime as a parent image
## final Dockerfile to build docker image
## version: 1.2
## Date: 18/10/2024
## last update: 25/10/2024
## last update: 29/10/2024
## last update: 20/11/2024 and workd
## last worked: 09/01/2025
## last worked: 04/02/2025
## worked with skip proccess, if faild any proccess, with security till VEP



##Download
## LoFtool_scores.txt
## environmenttmb.yml
## p2manta.yml
## vusprize1.yml

FROM ubuntu:20.04

# Set the working directory in the container
WORKDIR /usr/src/app
ENV DEBIAN_FRONTEND=noninteractive

#INSTALL dependencies
RUN apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get update && apt-get install -y \
        wget \
        python3 \
        python3-pip \
        bzip2 \
        gcc \
        make \
        unzip \
        curl \
	awscli \
        tzdata \
        git \
        bcftools \
        nano \
        less \
        autoconf \
        build-essential \
        cmake \
        g++ \
        gfortran \
        libcurl4-gnutls-dev \
        hdf5-tools \
        libboost-date-time-dev \
        libboost-program-options-dev \
        libboost-system-dev \
        libboost-filesystem-dev \
        libboost-iostreams-dev \
        libbz2-dev \
        libdeflate-dev \
        libhdf5-dev \
        libncurses-dev \
        liblzma-dev \
        pkg-config \
        zlib1g-dev \
        perl \
        cpanminus \
        libmysqlclient-dev \
        libxml2 \
        libexpat1-dev \
        libbz2-dev \
        liblzma-dev \
        libxml2-dev \
        zlib1g-dev  \
        && apt-get clean \
        && rm -rf /var/lib/apt/lists/*

# install delly
RUN cd /opt \
    && git clone --recursive https://github.com/dellytools/delly.git \
    && cd /opt/delly/ \
    && make STATIC=1 all \
    && make install

RUN pip install numpy tensorflow einops local_attention==1.0.0 pandas
# Update the package list
#RUN apt-get update
# Install OpenJDK 17
RUN apt-get update && apt-get install -y openjdk-17-jre-headless && apt-get clean
# Set JAVA_HOME environment variable
ENV JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64
ENV PATH="$JAVA_HOME/bin:${PATH}"

## install miniconda
RUN mkdir -p /opt/miniconda3 && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt/miniconda3/miniconda.sh && \
    bash /opt/miniconda3/miniconda.sh -b -u -p /opt/miniconda3 && \
    rm /opt/miniconda3/miniconda.sh
# Add conda to PATH
ENV PATH="/opt/miniconda3/bin:${PATH}" 
# Verify conda installation
RUN conda --version

# Install perl module
RUN apt-get update && \
    apt-get -y install gcc perl cpanminus libmysqlclient-dev libxml2 libexpat1-dev libbz2-dev liblzma-dev libxml2-dev zlib1g-dev && \
    apt-get clean

# Install the required Perl modules before installing BioPerl
RUN cpanm XML::Parser::PerlSAX && \
    cpanm XML::Twig && \
    cpanm LWP::UserAgent && \
    cpanm XML::DOM && \
    cpanm --force XML::LibXML && \
    cpanm --force XML::LibXML::Reader && \
    cpanm App::cpanminus Config::Tiny && \
    cpanm Module::Build && \
    cpanm Test::Warnings && \
    cpanm DBI && \
    cpanm DBD::mysql && \
    cpanm --force BioPerl && \
    cpanm --force Bio::Perl && \
    apt-get clean

# Clone and install Ensembl API
RUN rm -rf /usr/src/app/ensembl-vep && \ 
    git clone https://github.com/ensembl/ensembl-vep /usr/src/app/ensembl-vep && \
    cd /usr/src/app/ensembl-vep && \
    cpanm --installdeps . && \
    cpanm Module::Build && \
    perl INSTALL.pl && \
    perl INSTALL.pl -a p --PLUGINS all --PLUGINSDIR ./Plugins/

#Install fastQC from source
RUN wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip -O /tmp/fastqc.zip && \
    unzip /tmp/fastqc.zip -d /opt/ && \
    chmod +x /opt/FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc && \
    rm /tmp/fastqc.zip

#Install Fastp
RUN wget http://opengene.org/fastp/fastp.0.23.4 && \
    mv fastp.0.23.4 fastp && \
    chmod a+x ./fastp && \
    mv ./fastp /usr/local/bin

# Install MultiQC
RUN pip3 install multiqc

# Install BWA, SAMtools, Picard, and GATK using system packages or binaries
RUN apt-get install -y bwa samtools tabix bcftools
RUN wget https://github.com/broadinstitute/picard/releases/download/2.27.1/picard.jar -O picard.jar
RUN chmod +x picard.jar
   ## cp picard.jar /usr/local/bin

# Install GATK 4.5.00
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip -O gatk4500.zip && \
    unzip gatk4500.zip

COPY LoFtool_scores.txt /usr/src/app/ensembl-vep/Plugins/


# Install msisensor2
RUN git clone https://github.com/niu-lab/msisensor2.git && \
    cd msisensor2 && \
    chmod +x msisensor2

# Install pyTMB 
RUN git clone https://github.com/bioinfo-pf-curie/TMB.git
COPY environmenttmb.yml /tmp/
RUN conda env create -f /tmp/environmenttmb.yml

# Install MANTA
COPY p2manta.yml /tmp/
RUN conda env create -f /tmp/p2manta.yml

# Copy Databses for VEP into container
RUN cd /usr/src/app/ensembl-vep && \
    mkdir plugins_db

## install Mutserve
RUN mkdir mutserveTool && \
    cd mutserveTool && \
    curl -sL mutserve.vercel.app  | bash

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && mv nextflow /usr/local/bin/ \
    && apt-get update && apt-get install -y gpg graphviz

RUN pip install pysam pandas intervaltree
RUN pip install NanoPlot nanofilt
RUN git clone https://github.com/lh3/minimap2 && \
    cd minimap2 && \
    make
RUN git clone https://github.com/WGLab/NanoCaller.git && \
    conda env create -f NanoCaller/environment.yml && \
    chmod +x NanoCaller/NanoCaller
# Install MANTIS for Tumor - Normal MSI detection
RUN git clone https://github.com/OSU-SRLab/MANTIS.git

RUN conda create -n p3919 python=3.9.19 && \
    /bin/bash -c "source activate p3919 && pip install numpy pysam"
# Install CONTRA for CNV tumor only and Tumor-Normal
COPY CONTRA.v2.0.8.tar.gz /usr/src/app && \
    tar -xvzf CONTRA.v2.0.8.tar.gz && \
    rm CONTRA.v2.0.8.tar.gz
RUN apt-get install -y bedtools

RUN git clone https://github.com/mskcc/vcf2maf.git
RUN git clone https://github.com/oncokb/oncokb-annotator.git
RUN pip install xlsxwriter
COPY vusprize1.yml /tmp/
RUN conda env create -f /tmp/vusprize1.yml
RUN git clone https://github.com/danielhmahecha/VusPrize.git
WORKDIR /usr/src/app
# Copy the Nextflow script into the Docker image
RUN mkdir Validation_script
RUN apt-get update && apt-get install -y jq


#COPY icgeb15012025.nf /usr/src/app # for original validation
COPY VgenX.nf /usr/src/app

#COPY nextflow.config /usr/src/app
# Make port 80 available to the world outside this container
EXPOSE 80
# Define environment variable
ENV NAME rgenx
# Add Delly to PATH
#ENV PATH="/opt/delly/bin:${PATH}"
# Set entrypoint to use Snakemake
#ENTRYPOINT ["snakemake"]

# Default command that runs if no other command is specified
#CMD ["all", "--cores", "5"]
#ENV GPG_PASSPHRASE=10xl8$nhuUms340#rby7
# Default command: decrypt, run, re-encrypt, and clean up
CMD bash -c "\
#gpg --batch --yes --passphrase '$GPG_PASSPHRASE' -o main_securing3rd.nf -d main_securing3rd.nf.gpg && \
nextflow -log /usr/src/app/output/pipeline.log run VgenX.nf -c nextflow.config -params-file config.json -with-report /usr/src/app/output/pipeline_report -with-dag /usr/src/app/output/pipeline_DAG.png -with-trace /usr/src/app/output/pipeline_trace.txt -with-timeline /usr/src/app/output/pipeline_timeline.html && \
chmod -R 777 /usr/src/app/output"
#nextflow -log /usr/src/app/output/pipeline.log run VgenX.nf -c nextflow.config -params-file config.json -work-dir /usr/src/app/work -with-report /usr/src/app/output/pipeline_report -with-dag /usr/src/app/output/pipeline_DAG.png -with-trace /usr/src/app/output/pipeline_trace.txt -with-timeline /usr/src/app/output/pipeline_timeline.html && \

#gpg --batch --yes --passphrase '$GPG_PASSPHRASE' --symmetric --cipher-algo AES256 -o main_securing3rd.nf.gpg main_securing3rd.nf && \
#rm main_securing3rd.nf"
