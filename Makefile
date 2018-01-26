none:

install:
	curl -fsSL get.nextflow.io | bash

clean-logs:
	rm -f .nextflow.log.*

clean-output: 
	[ -d output ] && mv output oldoutput && rm -rf oldoutput &

clean-work: 
	[ -d work ] && mv work oldwork && rm -rf oldwork &

clean: clean-logs clean-work

clean-all: clean clean-output
