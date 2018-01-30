none:

# NextFlow setup & run commands
install:
	curl -fsSL get.nextflow.io | bash

run:
	./nextflow run pipeline.nf

clean-logs:
	rm -f .nextflow.log.*

clean-output:
	[ -d output ] && mv output oldoutput && rm -rf oldoutput &

clean-work:
	[ -d work ] && mv work oldwork && rm -rf oldwork &

clean: clean-logs clean-work

clean-all: clean clean-output


# other commands
docker:
	docker build -t stevekm/nextflow-demo .

# --privileged
docker-test:
	docker run --rm -ti stevekm/nextflow-demo
