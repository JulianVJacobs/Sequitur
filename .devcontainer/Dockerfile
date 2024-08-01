FROM alpine

RUN apk fix && apk --no-cache --update add git \
                                           build-base \
                                           python3-dev \
                                           py3-pip \
                                           linux-headers \
                                           perl \
                                           cmake \
                                           zlib-dev

RUN git clone https://github.com/ncbi/ncbi-vdb.git && \
    cd ncbi-vdb && ./configure && make
RUN git clone https://github.com/ncbi/sra-tools.git && \
    cd sra-tools && ./configure && make
RUN rm -r ncbi-vdb && rm -r sra-tools

ENV PATH="/root/ncbi-outdir/sra-tools/linux/gcc/x86_64/rel/bin/:${PATH}"

EXPOSE 41777