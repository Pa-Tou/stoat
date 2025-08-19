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

WORKDIR /stoat

COPY . /stoat

RUN git submodule update --init --recursive && \
    mkdir build \
    && cd build \
    && cmake .. \
    && make -j$(nproc)

ENV PATH=$PATH:/stoat/bin/

WORKDIR /home
