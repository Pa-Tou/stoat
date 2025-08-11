FROM ubuntu:22.04

# Éviter les interactions pendant l'installation
ENV DEBIAN_FRONTEND=noninteractive

# Mettre à jour le système et installer les dépendances de base
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    pkg-config \
    wget \
    libhts-dev \
    libboost-all-dev \
    libjansson-dev \
    protobuf-compiler \
    libprotoc-dev \
    libprotobuf-dev \
    valgrind \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /bin

RUN wget https://github.com/vgteam/vg/releases/download/v1.67.0/vg \
    && chmod +x vg

ENV PATH=$PATH:/bin/


WORKDIR /home

# Make sure that this gets rerun if the repo is updated
# From https://stackoverflow.com/questions/36996046/how-to-prevent-dockerfile-caching-git-clone
ADD https://api.github.com/repos/Plogeur/STOAT/git/refs/heads/main version.json

# Clone the STOAT C++ repository and set it as the working directory
RUN git clone --recursive https://github.com/Plogeur/STOAT \
    && cd STOAT \ 
    && git checkout stoat_cxx \
    && git submodule update --init --recursive \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make -j$(nproc)

ENV PATH=$PATH:/home/STOAT/bin/
