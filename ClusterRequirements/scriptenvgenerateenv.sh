#Start in a folder that contains only Infrared-master.tar.gz
module load cmake
module load gcc/11.3.0
module load gengetopt
module load python/3.10.2
module load python-build-bundle/2022a

git clone git@github.com:ViennaRNA/ViennaRNA.git
python -m venv ENV
source ENV/bin/activate
python -m pip install swig

cd ViennaRNA
git fetch
git checkout user-contrib
cd src
tar -xzf libsvm-3.25.tar.gz
tar -xjf dlib-19.23.tar.bz2
cd ..
autoreconf -fi
export PYTHON3=`which python3`; ./configure --without-python2 --without-perl --prefix=`pwd`
make
make install


cd ..
python -m pip install treedecomp
tar -xzvf Infrared-master.tar.gz
cd Infrared-master
python -m pip install .
cd ..
python -m pip install --no-index --upgrade pip
python -m pip install varnaapi
python -m pip install --no-index networkx
python -m pip install --no-index graphviz
python -m pip install pygraphviz
python -m pip install func_timeout
