none:

base:
	cd NGS580-base && docker build -t stevekm/ngs580-base .

hg19ref: base
	cd NGS580-hg19ref && docker build -t stevekm/ngs580-hg19ref .

main: hg19ref
	cd NGS580-main && docker build -t stevekm/ngs580-main .

build-test:
	cd test && docker build -t stevekm/test .

build: main

test:
	docker run --rm -ti stevekm/ngs580-hg19ref bash

test-test: build-test
	docker run --rm -ti stevekm/test bash
