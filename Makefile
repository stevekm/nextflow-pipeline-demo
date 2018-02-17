SHELL:=/bin/bash
none:

# ~~~~~ NEXTFLOW PIPELINE ~~~~~ #
# NextFlow setup & run commands
./nextflow:
	curl -fsSL get.nextflow.io | bash

install: ./nextflow

setup: 
	./generate-samplesheets.py example-data/ && ./update-samplesheets.py

# preprocessing pipeline
pre:  install
	./nextflow run preprocessing.nf -with-report "nextflow-preprocessing.html" -with-trace -with-timeline "timeline-preprocessing.html" -with-dag flowchart-preprocessing.png && \
	[ -f "trace.txt" ] && /bin/mv trace.txt trace-preprocessing.txt || :

prer: install
	./nextflow run preprocessing.nf -resume -with-report "nextflow-preprocessing.html" -with-trace -with-timeline "timeline-preprocessing.html" -with-dag flowchart-preprocessing.png && \
	[ -f "trace.txt" ] && /bin/mv trace.txt trace-preprocessing.txt || :

# exomes pipeline
ex: install
	./nextflow run exome.nf -with-report "nextflow-exome.html" -with-trace -with-timeline "timeline-exome.html" -with-dag flowchart-exome.png && \
	[ -f "trace.txt" ] && /bin/mv trace.txt trace-exome.txt || :

exr: install
	./nextflow run exome.nf -resume -with-report "nextflow-exome.html" -with-trace -with-timeline "timeline-exome.html" -with-dag flowchart-exome.png && \
	[ -f "trace.txt" ] && /bin/mv trace.txt trace-exome.txt || :

# exome-pairs pipeline
exp: install
	./nextflow run exome-pairs.nf -with-report "nextflow-exome-pairs.html" -with-trace -with-timeline "timeline-exome-pairs.html" -with-dag flowchart-exome-pairs.png && \
	[ -f "trace.txt" ] && /bin/mv trace.txt trace-exome-pairs.txt || :

expr: install
	./nextflow run exome-pairs.nf -resume -with-report "nextflow-exome-pairs.html" -with-trace -with-timeline "timeline-exome-pairs.html" -with-dag flowchart-exome-pairs.png && \
	[ -f "trace.txt" ] && /bin/mv trace.txt trace-exome-pairs.txt || :

# entire pipeline
ex-all: pre ex exp


#
# run: install
# 	./nextflow run pipeline.nf
#
# resume: install
# 	./nextflow run pipeline.nf -resume
#
# annot: install
# 	./nextflow run annotate.nf
#
#
# samples.fastq-raw.csv:
# 	./gather-fastqs.pl example-data/
#
# wes: install samples.fastq-raw.csv
# 	./nextflow run wes.nf -with-report "nextflow-wes.html" -with-trace -with-timeline "timeline-wes.html" -with-dag flowchart-wes.png && \
# 	[ -f "trace.txt" ] && /bin/mv trace.txt trace-wes.txt || :
#
# wes-r: install samples.fastq-raw.csv
# 	./nextflow run wes.nf -with-report "nextflow-wes.html" -with-trace -with-timeline "timeline-wes.html" -with-dag flowchart-wes.png -resume && \
# 	[ -f "trace.txt" ] && /bin/mv trace.txt trace-wes.txt || :
#
# wes-pairs:
# 	./nextflow run wes.nf -with-report "nextflow-wes-pairs.html" -with-trace -with-timeline "timeline-wes-pairs.html" -with-dag flowchart-wes-pairs.png && \
# 	[ -f "trace.txt" ] && /bin/mv trace.txt trace-wes-pairs.txt || :
#
# wes-all: wes wes-pairs


# ~~~~~ CLEANUP ~~~~~ #
clean-traces:
	rm -f trace.txt.* trace*.txt

clean-logs:
	rm -f .nextflow.log.*

clean-timelines:
	rm -f *.html

clean-reports:
	rm -f *.html

clean-flowcharts:
	rm -f *.png

clean-output:
	[ -d output ] && mv output oldoutput && rm -rf oldoutput &

clean-work:
	[ -d work ] && mv work oldwork && rm -rf oldwork &

clean: clean-logs clean-traces clean-timelines clean-reports clean-flowcharts

clean-all: clean clean-output clean-work

clean-wes: clean clean-work
	[ -d wes_output ] && mv wes_output wes_outputold && rm -rf wes_outputold &
	[ -d output-exomes ] && mv output-exomes output-exomesold && rm -rf output-exomesold &



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
