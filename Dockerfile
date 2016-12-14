FROM ubuntu:14.04

RUN mkdir -p /app
WORKDIR /app

RUN apt-get update
RUN apt-get install -y make rsync wget
RUN apt-get install -y git g++ libboost-all-dev libbz2-dev doxygen xsltproc docbook docbook-xsl docbook-xml autoconf automake autotools-dev
RUN mkdir -p /deps

# Install libzeep
RUN git clone https://github.com/mhekkel/libzeep.git /deps/libzeep ;\
    cd /deps/libzeep ;\
    git checkout tags/v3.0.3
# XXX: Workaround due to bug in libzeep's makefile
RUN sed -i '71s/.*/\t\$\(CXX\) \-shared \-o \$@ \-Wl,\-soname=\$\(SO_NAME\) \$\(OBJECTS\) \$\(LDFLAGS\)/' /deps/libzeep/makefile
WORKDIR /deps/libzeep
# XXX: Run ldconfig manually to work around a bug in libzeep's makefile
RUN make ; make install ; ldconfig

WORKDIR /app

COPY . /app

RUN ./autogen.sh && ./configure && make
