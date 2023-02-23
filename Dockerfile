# Minimal Docker image for ViralPrimerDesign using Ubuntu base
FROM ubuntu:20.04
MAINTAINER Niema Moshiri <niemamoshiri@gmail.com>

# set up environment + dependencies
RUN apt-get -qq update && apt-get -qq -y upgrade && \
    apt-get -qq install -y g++ python3 python3-pip wget && \

    # install Python packages
    pip3 install --no-cache-dir 'seaborn' && \

    # install MAFFT v7.505
    wget -qO- "https://mafft.cbrc.jp/alignment/software/mafft-7.505-without-extensions-src.tgz" | tar -zx && \
    cd mafft-*/core && \
    make clean && \
    make && \
    make install && \
    cd ../.. && \
    rm -rf mafft-* && \

    # install Primer3 v2.6.1
    wget -qO- "https://github.com/primer3-org/primer3/archive/refs/tags/v2.6.1.tar.gz" | tar -zx && \
    cd primer3-*/src && \
    make && \
    make install && \
    cd ../.. && \
    rm -rf primer3-* && \

    # install BLAST+ v2.13.0
    wget -qO- "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/ncbi-blast-2.13.0+-x64-linux.tar.gz" | tar -zx && \
    mv ncbi-blast-*/bin/* /usr/local/bin/ && \
    rm -rf ncbi-blast-* && \

    # clean up
    rm -rf /tmp/*
