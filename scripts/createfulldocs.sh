#! /bin/bash

# this script creates the full documentation, including python
# For this to work, gitlab will ask for credentials.
git clone https://git.physik.uni-wuerzburg.de/marqov/pymarqov.git
./scripts/compilefulldocs.sh
