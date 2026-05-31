FROM continuumio/miniconda3:4.12.0
WORKDIR /usr/src/app
RUN apt-get update && apt-get install -y \
    git \
    wget \
    curl \
    unzip \
    bzip2 \
    build-essential \
    autoconf \
    cmake \
    g++ \
    gfortran \
    libcurl4-gnutls-dev \
    libbz2-dev \
    libdeflate-dev \
    libhdf5-dev \
    libncurses-dev \
    liblzma-dev \
    pkg-config \
    zlib1g-dev \
    libxml2-dev \
    libexpat1-dev \
    libmysqlclient-dev \
    awscli \
    jq \
    nano \
    less \
    tzdata \
    graphviz \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
RUN apt-get update && apt-get install -y openjdk-17-jre-headless && apt-get clean
ENV JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64
ENV PATH="$JAVA_HOME/bin:${PATH}"
RUN conda env create -f /env/environment.yml 
RUN conda env create -f /env/environmenttmb.yml
RUN conda env create -f /env/p2manta.yml
RUN conda env create -f /env/vusprize1.yml
ENV PATH="/opt/conda/envs/wessom/bin:$PATH"
RUN curl -s https://get.nextflow.io | bash && mv nextflow /usr/local/bin/
COPY CONTRA.v2.0.8.tar.gz ./
RUN tar -xvzf CONTRA.v2.0.8.tar.gz && rm CONTRA.v2.0.8.tar.gz
RUN git clone https://github.com/mskcc/vcf2maf.git \
    && git clone https://github.com/oncokb/oncokb-annotator.git \
    && git clone https://github.com/WGLab/NanoCaller.git \
    && git clone https://github.com/OSU-SRLab/MANTIS.git \
    && git clone https://github.com/niu-lab/msisensor2.git && cd msisensor2 && chmod +x msisensor2 \
    && git clone https://github.com/danielhmahecha/VusPrize.git \
    && chmod +x NanoCaller/NanoCaller
RUN mkdir -p /opt/conda/envs/wessom/share/ensembl-vep/Plugins \
    && cp LoFtool_scores.txt /opt/conda/envs/wessom/share/ensembl-vep/Plugins/

RUN mkdir mutserveTool && cd mutserveTool && curl -sL mutserve.vercel.app | bash
RUN mkdir -p Validation_script output
COPY main.nf ./
EXPOSE 80
ENV NAME wessom

CMD bash -c "\
    nextflow -log /usr/src/app/output/pipeline.log run main.nf -c nextflow.config -params-file config.json \
    -with-report /usr/src/app/output/pipeline_report \
    -with-dag /usr/src/app/output/pipeline_DAG.png \
    -with-trace /usr/src/app/output/pipeline_trace.txt \
    -with-timeline /usr/src/app/output/pipeline_timeline.html && \
    chmod -R 777 /usr/src/app/output"
