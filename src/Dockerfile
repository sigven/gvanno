###################################################
# Stage 1 - docker container to build ensembl-vep #
###################################################
FROM ubuntu:18.04 as builder

# Update aptitude and install some required packages
# a lot of them are required for Bio::DB::BigFile
RUN apt-get update && apt-get -y install \
    build-essential \
    git \
    libpng-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    perl \
    perl-base \
    unzip \
    wget && \
    rm -rf /var/lib/apt/lists/*

# Setup VEP environment
ENV OPT /opt/vep
ENV OPT_SRC $OPT/src
ENV HTSLIB_DIR $OPT_SRC/htslib
ENV BRANCH release/105

# Working directory
WORKDIR $OPT_SRC

# Clone/download repositories/libraries
RUN if [ "$BRANCH" = "master" ]; \
    then export BRANCH_OPT=""; \
    else export BRANCH_OPT="-b $BRANCH"; \
    fi && \
    # Get ensembl cpanfile in order to get the list of the required Perl libraries
    wget -q "https://raw.githubusercontent.com/Ensembl/ensembl/$BRANCH/cpanfile" -O "ensembl_cpanfile" && \
    # Clone ensembl-vep git repository
    git clone $BRANCH_OPT --depth 1 https://github.com/Ensembl/ensembl-vep.git && chmod u+x ensembl-vep/*.pl && \
    # Clone ensembl-variation git repository and compile C code
    git clone $BRANCH_OPT --depth 1 https://github.com/Ensembl/ensembl-variation.git && \
    mkdir var_c_code && \
    cp ensembl-variation/C_code/*.c ensembl-variation/C_code/Makefile var_c_code/ && \
    rm -rf ensembl-variation && \
    chmod u+x var_c_code/* && \
    # Clone bioperl-ext git repository - used by Haplosaurus
    git clone --depth 1 https://github.com/bioperl/bioperl-ext.git && \
    # Download ensembl-xs - it contains compiled versions of certain key subroutines used in VEP
    wget https://github.com/Ensembl/ensembl-xs/archive/2.3.2.zip -O ensembl-xs.zip && \
    unzip -q ensembl-xs.zip && mv ensembl-xs-2.3.2 ensembl-xs && rm -rf ensembl-xs.zip && \
    # Clone/Download other repositories: bioperl-live is needed so the cpanm dependencies installation from the ensembl-vep/cpanfile file takes less disk space
    ensembl-vep/travisci/get_dependencies.sh && \
    # Only keep the bioperl-live "Bio" library
    mv bioperl-live bioperl-live_bak && mkdir bioperl-live && mv bioperl-live_bak/Bio bioperl-live/ && rm -rf bioperl-live_bak && \
    ## A lot of cleanup on the imported libraries, in order to reduce the docker image ##
    rm -rf Bio-HTS/.??* Bio-HTS/Changes Bio-HTS/DISCLAIMER Bio-HTS/MANIFEST* Bio-HTS/README Bio-HTS/scripts Bio-HTS/t Bio-HTS/travisci \
           bioperl-ext/.??* bioperl-ext/Bio/SeqIO bioperl-ext/Bio/Tools bioperl-ext/Makefile.PL bioperl-ext/README* bioperl-ext/t bioperl-ext/examples \
           ensembl-vep/.??* ensembl-vep/docker \
           ensembl-xs/.??* ensembl-xs/TODO ensembl-xs/Changes ensembl-xs/INSTALL ensembl-xs/MANIFEST ensembl-xs/README ensembl-xs/t ensembl-xs/travisci \
           htslib/.??* htslib/INSTALL htslib/NEWS htslib/README* htslib/test && \
    # Only keep needed kent-335_base libraries for VEP - used by Bio::DB::BigFile (bigWig parsing)
    mv kent-335_base kent-335_base_bak && mkdir -p kent-335_base/src && \
    cp -R kent-335_base_bak/src/lib kent-335_base_bak/src/inc kent-335_base_bak/src/jkOwnLib kent-335_base/src/ && \
    cp kent-335_base_bak/src/*.sh kent-335_base/src/ && \
    rm -rf kent-335_base_bak

# Setup bioperl-ext
WORKDIR bioperl-ext/Bio/Ext/Align/
RUN perl -pi -e"s|(cd libs.+)CFLAGS=\\\'|\$1CFLAGS=\\\'-fPIC |" Makefile.PL

# Install htslib binaries (for 'bgzip' and 'tabix')
# htslib requires the packages 'zlib1g-dev', 'libbz2-dev' and 'liblzma-dev'
WORKDIR $HTSLIB_DIR
RUN make install && rm -f Makefile *.c

# Compile Variation LD C scripts
WORKDIR $OPT_SRC/var_c_code
RUN make && rm -f Makefile *.c


###################################################
# Stage 2 - docker container to build ensembl-vep #
###################################################
FROM ubuntu:18.04

# Update aptitude and install some required packages
# a lot of them are required for Bio::DB::BigFile
RUN apt-get update && apt-get -y install \
    build-essential \
    cpanminus \
    curl \
    libmysqlclient-dev \
    libpng-dev \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    locales \
    openssl \
    perl \
    perl-base \
    unzip \
    wget \
    vim && \
    apt-get -y purge manpages-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Setup VEP environment
ENV OPT /opt/vep
ENV OPT_SRC $OPT/src
ENV PERL5LIB_TMP $PERL5LIB:$OPT_SRC/ensembl-vep:$OPT_SRC/ensembl-vep/modules
ENV PERL5LIB $PERL5LIB_TMP:$OPT_SRC/bioperl-live
ENV KENT_SRC $OPT/src/kent-335_base/src
ENV HTSLIB_DIR $OPT_SRC/htslib
ENV MACHTYPE x86_64
ENV DEPS $OPT_SRC
ENV PATH $OPT_SRC/ensembl-vep:$OPT_SRC/var_c_code:$PATH
ENV LANG_VAR en_US.UTF-8

# Create vep user
RUN useradd -r -m -U -d "$OPT" -s /bin/bash -c "VEP User" -p '' vep && usermod -a -G sudo vep && mkdir -p $OPT_SRC
USER vep

# Copy downloaded libraries (stage 1) to this image (stage 2)
COPY --chown=vep:vep --from=builder $OPT_SRC $OPT_SRC
#############################################################

# Change user to root for the following complilations/installations
USER root

# Install bioperl-ext, faster alignments for haplo (XS-based BioPerl extensions to C libraries)
WORKDIR $OPT_SRC/bioperl-ext/Bio/Ext/Align/
RUN perl Makefile.PL && make && make install && rm -f Makefile*

# Install ensembl-xs, faster run using re-implementation in C of some of the Perl subroutines
WORKDIR $OPT_SRC/ensembl-xs
RUN perl Makefile.PL && make && make install && rm -f Makefile* cpanfile

WORKDIR $OPT_SRC
# Install/compile more libraries
RUN ensembl-vep/travisci/build_c.sh && \
    # Remove unused Bio-DB-HTS files
    rm -rf Bio-HTS/cpanfile Bio-HTS/Build.PL Bio-HTS/Build Bio-HTS/_build Bio-HTS/INSTALL.pl && \
    # Install ensembl perl dependencies (cpanm)
    cpanm --installdeps --with-recommends --notest --cpanfile ensembl_cpanfile . && \
    cpanm --installdeps --with-recommends --notest --cpanfile ensembl-vep/cpanfile . && \
    # Delete bioperl and cpanfiles after the cpanm installs as bioperl will be reinstalled by the INSTALL.pl script
    rm -rf bioperl-live ensembl_cpanfile ensembl-vep/cpanfile && \
    # Configure "locale", see https://github.com/rocker-org/rocker/issues/19
    echo "$LANG_VAR UTF-8" >> /etc/locale.gen && locale-gen en_US.utf8 && \
    /usr/sbin/update-locale LANG=$LANG_VAR && \
    # Copy htslib executables. It also requires the packages 'zlib1g-dev', 'libbz2-dev' and 'liblzma-dev'
    cp $HTSLIB_DIR/bgzip $HTSLIB_DIR/tabix $HTSLIB_DIR/htsfile /usr/local/bin/

ENV LC_ALL $LANG_VAR
ENV LANG $LANG_VAR

# Switch back to vep user
USER vep
ENV PERL5LIB $PERL5LIB_TMP


USER vep
ENV PERL5LIB $PERL5LIB_TMP

WORKDIR /
ADD loftee_1.0.3.tgz $OPT/src/ensembl-vep/modules
ADD UTRannotator.tgz $OPT/src/ensembl-vep/modules
RUN wget -q "https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/105/NearestExonJB.pm" -O $OPT/src/ensembl-vep/modules/NearestExonJB.pm


# Final steps
WORKDIR $OPT_SRC/ensembl-vep
# Update bash profile
RUN echo >> $OPT/.profile && \
    echo PATH=$PATH:\$PATH >> $OPT/.profile && \
    echo export PATH >> $OPT/.profile && \
    # Run INSTALL.pl and remove the ensemb-vep tests and travis
    ./INSTALL.pl -a a -l -n && rm -rf t travisci .travis.yml

######################################################
# Stage 3 - adding dependencies for gvanno analysis  #
######################################################

USER root

RUN apt-get update \
  && apt-get install -y python3-pip python3-dev \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 install --upgrade pip

RUN apt-get update && apt-get -y install apache2 apt-utils build-essential cpanminus curl git libmysqlclient-dev libpng-dev libssl-dev manpages mysql-client openssl perl perl-base unzip vim wget sudo
# install ensembl dependencies
RUN cpanm Test::Object PPI::Document Task::Weaken Test::SubCalls Test::Object DBI DBD::mysql Archive::Zip Perl::Critic Set::IntervalTree
RUN apt-get update && apt-get install apt-transport-https

USER root
WORKDIR /

ENV PACKAGE_BIO="libhts2 bedtools"
ENV PACKAGE_DEV="gfortran gcc-multilib autoconf liblzma-dev libncurses5-dev libblas-dev liblapack-dev libssh2-1-dev libxml2-dev vim libssl-dev libcairo2-dev libbz2-dev libcurl4-openssl-dev"
#ENV PYTHON_MODULES="numpy==1.19.2 cython==0.29.21 scipy==1.5.3 pandas==1.1.3 cyvcf2==0.20.9 toml==0.10.1"
ENV PYTHON_MODULES="numpy==1.19.2 cython==0.29.21 scipy==1.5.3 pandas==1.1.3 cyvcf2==0.20.9"
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		nano ed locales vim-tiny fonts-texgyre \
    $PACKAGE_DEV $PACKAGE_BIO \
	&& rm -rf /var/lib/apt/lists/*

RUN apt-get autoremove

## Install vcfanno version v0.3.2
RUN wget https://github.com/brentp/vcfanno/releases/download/v0.3.2/vcfanno_linux64 && \
    mv vcfanno_linux64 vcfanno && \
    mv vcfanno /usr/local/bin && \
    chmod 755 /usr/local/bin/vcfanno

## Install Ensembl's Vcf-validator
RUN wget https://github.com/EBIvariation/vcf-validator/releases/download/v0.9.3/vcf_validator_linux && \
mv vcf_validator_linux /usr/local/bin/vcf_validator && \
chmod 755 /usr/local/bin/vcf_validator

USER root
WORKDIR /

RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
RUN bunzip2 -dc samtools-1.10.tar.bz2 | tar xvf -
RUN cd samtools-1.10 && ./configure --prefix=/usr/local/bin && make -j && make install

WORKDIR /

## Install tools used for compilation
RUN sudo -H pip install --upgrade pip
RUN sudo -H pip install -U setuptools
RUN sudo -H pip install $PYTHON_MODULES

RUN apt-get update \
  && apt-get install -y --no-install-recommends libpq-dev libxt-dev libudunits2-dev

USER root
WORKDIR /

## vt - variant tool set - use conda version
## primary use in GVANNO/PCGR/CPSR: decomposition of multiallelic variants in a VCF file
RUN wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh \
    && chmod 0755 miniconda.sh
RUN ["/bin/bash", "-c", "/miniconda.sh -b -p /conda"]
RUN rm miniconda.sh

# update conda & install vt
RUN /conda/bin/conda update conda
RUN /conda/bin/conda install -c bioconda vt

## Clean Up
RUN apt-get clean autoclean
RUN rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN rm -rf /var/lib/{dpkg,cache,log}

VOLUME /workdir
WORKDIR /workdir/
USER root
RUN mkdir /data && chmod 777 /data
WORKDIR /data
VOLUME /data

USER root
WORKDIR /
#RUN rm -f nlopt-2.4.2.tar.gz
RUN rm -rf $HOME/src/ensembl-vep/t/
RUN rm -f $HOME/src/v335_base.tar.gz
RUN rm -f $HOME/src/release-1-6-924.zip
RUN rm -rf /samtools-1.10.tar.bz2
RUN rm -f /conda/bin/python

ADD gvanno.tgz /
ENV PATH=$PATH:/conda/bin:/gvanno
ENV PYTHONPATH=:/gvanno/lib:${PYTHONPATH}