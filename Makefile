main:
	@echo "TODO:probably change the target here"
	python FuzzTree.py #>> output.txt

dependencies:
	conda install -c conda-forge -c bioconda 'infrared' viennarna jupyter matplotlib pip networkx graphviz pygraphviz
	pip install varnaapi

tests:
	@echo "TODO"

clean:
	rm *.png
	rm *.pdf
	rm output.txt
