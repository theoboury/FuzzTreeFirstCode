main:
	@echo "TODO:probably change the target here"
	python TestFuzzTree.py >> output.txt

dependencies:
	conda install -c conda-forge -c bioconda 'infrared' viennarna jupyter matplotlib pip networkx graphviz pygraphviz
	pip install varnaapi
	pip install func_timeout

tests:
	python TestFuzzTree.py
clean:
	rm -f *.png
	rm -R -f __pycache__
	rm -f *.pdf
	rm -f output.txt
