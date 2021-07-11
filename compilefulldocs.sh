#! /bin/bash

# this script compiles the full documentation, including python.
# It assumes that the sources are already present at the specified positions.
rm pymarqov/doc/conf.py
cd doc
cp -rn ../pymarqov/doc python
echo "sys.path.insert(0, os.path.abspath('../pymarqov/'))" >> conf.py
doxygen Doxyfile-marqov-breathe.cfg
make html
