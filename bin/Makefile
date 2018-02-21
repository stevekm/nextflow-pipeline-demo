none:



venv/bin/activate:
	export PYTHONPATH= && \
	# unset PYTHONPATH && \
	virtualenv venv --no-site-packages &&\
	ln -fs venv/bin/activate

setup: venv/bin/activate
	export PYTHONPATH= && \
	source venv/bin/activate && \
	pip install -r requirements.txt


clean:
	rm -rf venv
	unlink activate