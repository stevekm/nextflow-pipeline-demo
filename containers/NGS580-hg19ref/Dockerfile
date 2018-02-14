FROM stevekm/ngs580-base

MAINTAINER Stephen M. Kelly

ENV HG19_GENOME_FA_MD5 c1ddcc5db31b657d167bea6d9ff354f9
ENV HG19_FA /opt/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
RUN cd opt/ && wget -r --no-parent -e robots=off -nH --cut-dirs=4 https://genome.med.nyu.edu/results/external/NYU/snuderllab/ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/
RUN find opt/ref/ -type f -name "index.html*" -delete
# RUN download_md5="$(md5sum /opt/iGenomes/Homo_sapiens2/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa | cut -f1 -d ' ')" && [ "${download_md5}" != "${HG19_GENOME_FA_MD5}" ] && echo "the md5sum doesnt match" && exit 1 || echo "the md5sums match"
# RUN mv /opt/iGenomes/Homo_sapiens2 /opt/iGenomes/Homo_sapiens
# RUN [ ! -f "${HG19_FA}" ] && echo "hg19 reference file not found at location: ${HG19_FA}" && exit 1 || echo "hg19 reference file is present"
