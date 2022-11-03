main:
	@echo "TODO:probably change the target here"
	python TestFuzzTree.py >> output.txt

dependencies:
	conda install -c conda-forge -c bioconda 'infrared' viennarna jupyter matplotlib pip networkx graphviz pygraphviz
	pip install varnaapi

tests:
	python TestFuzzTree.py
clean:
	rm *.png
	rm *.pdf
	rm output.txt
