SHELL:=/bin/bash
none:

# ~~~~~ NEXTFLOW PIPELINE ~~~~~ #
# NextFlow setup & run commands
./nextflow:
	curl -fsSL get.nextflow.io | bash

install: ./nextflow

run: install
	./nextflow run pipeline.nf

resume: install
	./nextflow run pipeline.nf -resume

annot: install
	./nextflow run annotate.nf


samples.fastq-raw.csv:
	./gather-fastqs.pl example-data/

wes: install samples.fastq-raw.csv
	./nextflow run wes.nf -with-report -with-trace -with-timeline -with-dag flowchart.png





# ~~~~~ CLEANUP ~~~~~ #
clean-traces:
	rm -f trace.txt.*

clean-logs:
	rm -f .nextflow.log.*

clean-timelines:
	rm -f timeline.html.*

clean-reports:
	rm -f report.html.*
	rm -f nextflow-report.html.*

clean-flowcharts:
	rm -f flowchart.png.*

clean-output:
	[ -d output ] && mv output oldoutput && rm -rf oldoutput &

clean-work:
	[ -d work ] && mv work oldwork && rm -rf oldwork &

clean: clean-logs clean-traces clean-timelines clean-reports clean-flowcharts

clean-all: clean clean-output clean-work

clean-wes: clean
	[ -d wes_output ] && mv wes_output wes_outputold && rm -rf wes_outputold &
	[ -d output-exomes ] && mv output-exomes output-exomesold && rm -rf output-exomesold &



# ~~~~~ SETUP ~~~~~ #
# DOCKER

# REFERENCE FILES
HG19_GENOME_FA_MD5:=c1ddcc5db31b657d167bea6d9ff354f9

ref:
	mkdir ref

ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa: ref
	wget -r --no-parent -e robots=off -nH --cut-dirs=4 https://genome.med.nyu.edu/results/external/NYU/snuderllab/ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/ && \
	bin/compare_md5.sh ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa "$(HG19_GENOME_FA_MD5)"


all: ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa

references: all
	find ref/ -type f -name "index.html*" -delete

clean-ref:
	[ -d ref ] && \
	mv ref oldref && \
	mkdir ref && \
	touch ref/.gitkeep \
	&& rm -rf oldref &
