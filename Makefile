main:
	@echo "TODO:probably change the target here"
	conda activate FuzzTreeEnv
	python TestFuzzTree.py >> output.txt
	conda deactivate

env:
	
	conda create --name FuzzTreeEnv

dependencies:
	conda activate FuzzTreeEnv
	conda install -c conda-forge -c bioconda 'infrared' viennarna jupyter matplotlib pip networkx graphviz pygraphviz
	pip install varnaapi
	pip install func_timeout
	conda deactivate
tests:
	conda activate FuzzTreeEnv
	python TestFuzzTree.py
	conda deactivate
clean:
	rm -f *.png
	rm -R -f __pycache__
	rm -f *.pdf
	rm -f output.txt
