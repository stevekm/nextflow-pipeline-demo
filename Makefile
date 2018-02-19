SHELL:=/bin/bash
none:

# ~~~~~ NEXTFLOW PIPELINE ~~~~~ #
# NextFlow setup & run commands
./nextflow:
	curl -fsSL get.nextflow.io | bash

install: ./nextflow

setup: 
	./generate-samplesheets.py example-data/ && ./update-samplesheets.py

# exomes pipeline
exome: install
	./nextflow run exome.nf -with-report "nextflow-exome.html" -with-trace -with-timeline "timeline-exome.html" -with-dag flowchart-exome.dot && \
	[ -f "trace.txt" ] && /bin/mv trace.txt trace-exome.txt || :

exomer: install
	./nextflow run exome.nf -resume -with-report "nextflow-exome.html" -with-trace -with-timeline "timeline-exome.html" -with-dag flowchart-exome.dot && \
	[ -f "trace.txt" ] && /bin/mv trace.txt trace-exome.txt || :

# entire pipeline
run: exome


# ~~~~~ CLEANUP ~~~~~ #
clean-traces:
	rm -f trace.txt.* trace*.txt

clean-logs:
	rm -f .nextflow.log.*

clean-timelines:
	rm -f *.html

clean-reports:
	rm -f *.html nextflow.html*

clean-flowcharts:
	rm -f *.png

clean-output:
	[ -d output ] && mv output oldoutput && rm -rf oldoutput &
	[ -d output-exome ] && mv output-exome oldoutput-exome && rm -rf oldoutput-exome &
	

clean-work:
	[ -d work ] && mv work oldwork && rm -rf oldwork &

clean: clean-logs clean-traces clean-timelines clean-reports clean-flowcharts

clean-all: clean clean-output clean-work 
	[ -d .nextflow ] && mv .nextflow .nextflowold && rm -rf .nextflowold &
	rm -f .nextflow.log

# ~~~~~ SETUP ~~~~~ #
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
