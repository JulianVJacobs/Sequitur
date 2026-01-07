FROM alpine

RUN apk fix && apk --no-cache --update add git \
                                           build-base \
                                           python3-dev \
                                           py3-pip \
                                           linux-headers \
                                           perl \
                                           cmake \
                                           zlib-dev \
                                           autoconf \
                                           automake \
                                           make \
                                           gcc \
                                           musl-dev \
                                           perl \
                                           bash \
                                           bzip2-dev \
                                           xz-dev \
                                           curl-dev \
                                           openssl-dev \
                                           ncurses-dev \
                                           libdeflate-dev \
                                           rustup \
                                           openssh-client

RUN rustup-init -y && \
    source $HOME/.cargo/env && \
    rustup component add rustfmt

# RUN git clone https://github.com/ncbi/ncbi-vdb.git && \
#     cd ncbi-vdb && ./configure && make
# RUN git clone https://github.com/ncbi/sra-tools.git && \
#     cd sra-tools && ./configure && make
# RUN rm -r ncbi-vdb && rm -r sra-tools

RUN wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2 && \
    tar -xvjf samtools-1.20.tar.bz2 && rm samtools-1.20.tar.bz2 && \
    cd samtools-1.20 && ./configure && make && make install

RUN git clone --depth 1 https://github.com/seqan/seqan.git /root/seqan && \
    mkdir -p /root/mason2 && cd /root/mason2 && \
    cmake /root/seqan -DCMAKE_BUILD_TYPE=Release -DSEQAN_BUILD_SYSTEM=APP:mason2 && \
    make -j"$(nproc)" && \
    cp bin/mason_* /usr/local/bin/ && \
    rm -rf /root/seqan /root/mason2

ENV SSH_DIR=/root/.ssh

EXPOSE 41777