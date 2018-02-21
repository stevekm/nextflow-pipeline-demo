# Installation script for Conpair program

none:

Conpair:
	git clone https://github.com/nygenome/Conpair.git

install: Conpair




Conpair/venv/bin/activate:
	module purge && module load local && module load python/2.7.3 && \
	export PYTHONPATH= && \
	# unset PYTHONPATH && \
	virtualenv Conpair/venv --no-site-packages &&\
	ln -fs Conpair/venv/bin/activate Conpair/activate

setup: install Conpair/venv/bin/activate
	module purge && module load local && module load python/2.7.3 && \
	export PYTHONPATH= && \
	source Conpair/venv/bin/activate && \
	pip install -r Conpair.requirements.txt




# 94dc877fda41269563607bf614786bc7  picard.jar
Conpair/picard.jar:
	wget -O Conpair/picard.jar https://github.com/broadinstitute/picard/releases/download/2.17.10/picard.jar

Conpair/data/genomes:
	mkdir -p Conpair/data/genomes

# 0ce84c872fc0072a885926823dcd0338  human_g1k_v37.fa
Conpair/data/genomes/human_g1k_v37.fa: Conpair/data/genomes
	-wget -qO- ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz | gunzip > Conpair/data/genomes/human_g1k_v37.fa

# 639401c27993c476e59c0dd0aaa1fd1f  human_g1k_v37.dict
Conpair/data/genomes/human_g1k_v37.dict: Conpair/data/genomes/human_g1k_v37.fa Conpair/picard.jar
	java -jar Conpair/picard.jar CreateSequenceDictionary  R=Conpair/data/genomes/human_g1k_v37.fa O=Conpair/data/genomes/human_g1k_v37.dict

# 772484cc07983aba1355c7fb50f176d4  human_g1k_v37.fa.fai
Conpair/data/genomes/human_g1k_v37.fa.fai: Conpair/data/genomes/human_g1k_v37.fa
	module purge && module load local && module load samtools/1.3 && \
	samtools faidx Conpair/data/genomes/human_g1k_v37.fa

ref: Conpair/data/genomes/human_g1k_v37.fa Conpair/data/genomes/human_g1k_v37.dict Conpair/data/genomes/human_g1k_v37.fa.fai




clean-ref:
	rm -f Conpair/human_g1k_v37.*
	rm -f Conpair/data/genomes/*

clean:
	rm -rf Conpair/venv &
	unlink Conpair/activate
