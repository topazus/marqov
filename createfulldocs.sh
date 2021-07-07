#! /bin/bash

# this script creates the full documentation, including python
# For this to work, gitlab will ask for credentials.
git clone https://git.physik.uni-wuerzburg.de/marqov/pymarqov.git
rm pymarqov/doc/conf.py
cd doc
cp -rn ../pymarqov/doc python
echo "sys.path.insert(0, os.path.abspath('../pymarqov/'))" >> conf.py
doxygen Doxyfile-marqov-breathe.cfg
make html
