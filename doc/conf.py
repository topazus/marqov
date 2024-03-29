# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('..'))


# -- Project information -----------------------------------------------------

project = 'pyMARQOV/MARQOV'
copyright = '2021-2022, The MARQOV project'
author = 'Manuel Schrauth, Florian Goth'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon', 'breathe', 'sphinxcontrib.plantuml', 'sphinx.ext.mathjax']

# from the CXX version
autosummary_generate = True

# For PlantUML
plantuml_output_format = 'svg'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_css_files = ['custom.css']
html_logo = "./_static/logo.png"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

mathjax_path='https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML'

# configure napoleon extension; makes Sphinx understand NumPy or Google style docstrings!
napoleon_google_docstring = False
napoleon_numpy_docstring = True

# set up breathe only if doxygen config is present.
if os.path.exists('./doxygen/xml'):
    print("Using doxygen!")
    # CXX
    breathe_projects = { "Marqov": "./doxygen/xml" }

    breathe_default_project = "Marqov"
