FROM alpine

RUN apk fix && apk --no-cache --update add git build-base python3-dev py3-pip linux-headers perl cmake

RUN git clone https://github.com/ncbi/ncbi-vdb.git ncbi-vdb && \
    cd ncbi-vdb && ./configure && make
RUN git clone https://github.com/ncbi/sra-tools.git sra-tools && \
    cd sra-tools && ./configure && make
RUN rm -r ncbi-vdb && rm -r sra-tools