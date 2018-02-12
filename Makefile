none:

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

clean-logs:
	rm -f .nextflow.log.*

clean-output:
	[ -d output ] && mv output oldoutput && rm -rf oldoutput &

clean-work:
	[ -d work ] && mv work oldwork && rm -rf oldwork &

clean: clean-logs clean-work

clean-all: clean clean-output


# other commands
build:
	docker build -t stevekm/nextflow-demo .

build-clean:
	docker build --no-cache -t stevekm/nextflow-demo .

# --privileged
build-test:
	docker run --rm -ti stevekm/nextflow-demo bash

docker-demo:
	docker run --privileged --rm -ti debian:jessie /bin/bash
