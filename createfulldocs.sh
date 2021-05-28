#! /bin/bash

# this script creates the full documentation, including python
# For this to work, gitlab wil lask for credentials.
git clone https://git.physik.uni-wuerzburg.de/marqov/jupyter.git
rm jupyter/pymarqov/doc/conf.py
cd doc
cp -rn ../jupyter/pymarqov/doc python
echo "sys.path.insert(0, os.path.abspath('../jupyter/pymarqov/'))" >> conf.py
doxygen Doxyfile-marqov-breathe.cfg
make html
