dependencies:
	conda activate FuzzTreeEnv
	conda install -c conda-forge -c bioconda 'infrared' viennarna jupyter matplotlib pip networkx graphviz pygraphviz
	pip install varnaapi
	pip install func_timeout
	conda deactivate
tests:
	python3 Launcher.py --task create_patterns_and_targets
	python Launcher.py --task launch_sanity_test
clean:
	rm -f *.png
	rm -R -f __pycache__
	rm -f *.pdf
	rm -f output.txt
