# demo container for pipeline dependencies
FROM centos:7

MAINTAINER Stephen M. Kelly

RUN yum update -y

RUN yum install -y wget
RUN yum install -y git
RUN yum install -y make
RUN yum install -y gcc
RUN yum install -y zlib
RUN yum install -y zlib-devel
RUN yum install -y bzip2
RUN yum install -y bzip2-libs
RUN yum install -y bzip2-devel
RUN yum install -y gcc-c++
RUN yum install -y unzip
RUN yum install -y xz-devel

RUN mkdir /opt/bin
ENV PATH="/opt/bin:${PATH}"

# ~~~~~ hstlib ~~~~~ #
RUN cd /opt/bin && wget https://github.com/samtools/htslib/releases/download/1.3.1/htslib-1.3.1.tar.bz2 && tar xvfj htslib-1.3.1.tar.bz2 && rm -f htslib-1.3.1.tar.bz2
RUN cd /opt/bin/htslib-1.3.1 && make

# ~~~~~ SAMTOOLS ~~~~~ #
RUN cd /opt/bin && wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && tar xvfj samtools-1.3.1.tar.bz2 && rm -f samtools-1.3.1.tar.bz2
RUN cd /opt/bin/samtools-1.3.1 && ./configure --without-curses && make && make install
# RUN cd /opt/bin && wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2 && tar xvfj samtools-1.7.tar.bz2 && rm -f samtools-1.7.tar.bz2
# RUN cd /opt/bin/samtools-1.7 && ./configure --without-curses && make && make install
# RUN exit 1

# # ~~~~~ vcflib ~~~~~ #
# RUN cd /opt/bin && git clone --recursive https://github.com/ekg/vcflib.git && cd vcflib && make

# ~~~~~ bcftools ~~~~~ #
RUN cd /opt/bin && wget https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2 && tar xvfj bcftools-1.3.1.tar.bz2 && cd bcftools-1.3.1 && make

# ~~~~~ IGV ~~~~~ #
RUN cd /opt/bin && wget http://data.broadinstitute.org/igv/projects/downloads/2.3/IGV_2.3.81.zip -O tmp && mv tmp IGV_2.3.81.zip && unzip IGV_2.3.81.zip && rm -f IGV_2.3.81.zip

# ~~~~~ PANDOC ~~~~~ #
RUN cd /opt/bin && wget https://s3.amazonaws.com/rstudio-buildtools/pandoc-1.13.1.zip && unzip pandoc-1.13.1.zip && rm -f pandoc-1.13.1.zip


# ~~~~~ hg19 reference fasta ~~~~~ #
# backup URL here (much larger) ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg19/Homo_sapiens_UCSC_hg19.tar.gz
ENV HG19_GENOME_FA_MD5 c1ddcc5db31b657d167bea6d9ff354f9
ENV HG19_FA /opt/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
RUN mkdir /opt/iGenomes && cd /opt/iGenomes && wget -q https://s3.amazonaws.com/reftransdata/Homo_sapiens_UCSC_hg19_small.tar.gz -O- | tar -zxvf - Homo_sapiens2/UCSC/hg19/Sequence/WholeGenomeFasta
RUN download_md5="$(md5sum /opt/iGenomes/Homo_sapiens2/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa | cut -f1 -d ' ')" && [ "${download_md5}" != "${HG19_GENOME_FA_MD5}" ] && echo "the md5sum doesnt match" && exit 1 || echo "the md5sums match"
RUN mv /opt/iGenomes/Homo_sapiens2 /opt/iGenomes/Homo_sapiens
RUN [ ! -f "${HG19_FA}" ] && echo "hg19 reference file not found at location: ${HG19_FA}" && exit 1 || echo "hg19 reference file is present"


# # ~~~~~ MSI SENSOR ~~~~~ #
ENV HG19_MICROSATELLITES /opt/bin/msisensor/hg19_microsatellites.list
ENV SAMTOOLS_ROOT /opt/bin/samtools-1.3.1
# ENV SAMTOOLS_ROOT /opt/bin/samtools-1.7
# RUN cd /opt/bin/samtools-1.3.1 && ln -s htslib-1.3.1/htslib
RUN cd /opt/bin && git clone https://github.com/ding-lab/msisensor.git
RUN cd /opt/bin/msisensor && make
# RUN cd /opt/bin/msisensor && ./msisensor scan -d "${HG19_FA}" -o "${HG19_MICROSATELLITES}"
# RUN [ ! -e "${HG19_MICROSATELLITES}" ] && echo "hg19 microsatellites list file is not found at location: ${HG19_MICROSATELLITES}" && exit 1
# RUN cd /opt/bin/msisensor/test && bash run.sh
