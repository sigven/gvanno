############################################################
# Dockerfile to build Germline Variant Annotator (gvanno)
# Main software components:
#  1. Variant Effect Predictor (VEP)
#  2. vcfanno
#  3. custom scripts (gvanno.tgz)
############################################################

# use Debian as base image
FROM debian:latest

# set non interactive
ENV DEBIAN_FRONTEND=noninteractive
ENV PACKAGE_UTILS="apt-utils wget curl nano unzip"
ENV PACKAGE_BIO="tabix samtools libhts1 bedtools"
ENV PACKAGE_DEV="perl debconf-utils build-essential gfortran python-dev python-pip gcc-multilib autoconf zlib1g-dev git libncurses5-dev libblas-dev liblapack-dev cpanminus"
ENV PYTHON_MODULES="numpy cython scipy bx-python pyvcf cyvcf biopython"

USER root
WORKDIR /

# Update the repository sources list
RUN apt-get update

# Install tools used for compilation
RUN apt-get install --yes --no-install-recommends $PACKAGE_UTILS $PACKAGE_DEV $PACKAGE_BIO
RUN pip install -U setuptools
RUN pip install $PYTHON_MODULES

# FROM VCFANNO DOCKERFILE: Install vcfanno version 0.0.11
RUN wget https://github.com/brentp/vcfanno/releases/download/v0.0.11/vcfanno_0.0.11_linux_amd64.tar.gz && \
    tar xvzf vcfanno_0.0.11_linux_amd64.tar.gz && \
    mv vcfanno_0.0.11_linux_amd64/vcfanno /usr/local/bin && \
    rm vcfanno_0.0.11_linux_amd64.tar.gz && \
    rm -rf vcfanno_0.0.11_linux_amd64


# Install VEP
# https://gist.github.com/ckandoth/57d189f018b448774704d3b2191720a6
RUN cpanm File::ShareDir::Install \
 && cpanm Data::UUID \
 && cpanm autodie \
 && cpanm Test::Fatal \
 && cpanm IPC::System::Simple \
 && cpanm File::Fetch \
 && cpanm File::Copy::Recursive \
 && cpanm LWP::Simple \
 && cpanm Bio::Root::Version \
 && cpanm Archive::Extract \
 && cpanm Archive::Zip \
 && cpanm CGI \
 && cpanm DBI \
 && cpanm Archive::Tar

# && cpanm Bio::DB::HTS \

 VOLUME        /usr/local/share/vep/data
 ENV VEP_DATA="/usr/local/share/vep/data"
 ENV VEP_DATA_DOCKER="/usr/local/share/vep/data"
 ENV VEP_PATH="/usr/local/share/vep/vep"
 ENV VEP_PLUGINS="/usr/local/share/vep/Plugins"
 ENV VEP_INSTDIR="/ensembl-tools/scripts/variant_effect_predictor/"
 ENV PERL5LIB=$PERL5LIB:$VEP_PATH:$VEP_PLUGINS
 ENV PATH=$PATH:$VEP_PATH/htslib
 ENV SPECIES="homo_sapiens"
 ENV ASSEMBLY="GRCh37"
 ENV VEPPLUGIN="LoF,TSSDistance"
 RUN git clone -b release/85 https://github.com/Ensembl/ensembl-tools.git
 WORKDIR $VEP_INSTDIR
 RUN printf 'y\n' | perl INSTALL.pl --AUTO a --SPECIES $SPECIES --ASSEMBLY $ASSEMBLY --PLUGINS $VEPPLUGIN --DESTDIR $VEP_PATH --CACHEDIR $VEP_DATA
 RUN cp *.pl /usr/local/bin/
 #RUN git clone -b master --depth=1 https://github.com/samtools/htslib.git
 #WORKDIR $VEP_INSTDIR/htslib
 #RUN make
 #WORKDIR $VEP_INSTDIR
 #RUN git clone https://github.com/Ensembl/Bio-HTS.git biodbhts
 ############### FIX FOR VERSION 84 ##############
 WORKDIR /


 ADD gvanno.tgz /
 ENV PATH=$PATH:/gvanno
 #ENV PATH=$PATH:/htslib
 ENV PYTHONPATH=:/gvanno/lib:${PYTHONPATH}
 ENV VCFANNO_DATA_DOCKER="/data"

# Clean Up
RUN apt-get clean autoclean
RUN apt-get autoremove -y
RUN rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN rm -rf /var/lib/{dpkg,cache,log}

VOLUME /workdir
WORKDIR /workdir/


### FROM VCFANNO DOCKERFILE
USER root
RUN mkdir /data && chmod 777 /data
WORKDIR /data
VOLUME /data
WORKDIR /
