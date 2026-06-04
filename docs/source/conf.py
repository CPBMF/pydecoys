# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import sys
from pathlib import Path

sys.path.insert(0, str(Path('..', '..').resolve()))


project = 'pydecoys'
copyright = \
    '2026, CPBMF and INCT-TB, PUCRS, Porto Alegre, Brazil, Bruno Maestri A Becker'

author = 'Bruno Maestri A Becker'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'sphinx_rtd_theme'
]

templates_path = ['_templates']
exclude_patterns = []

add_module_names = False

autodoc_typehints = "description"

napoleon_use_rtype = True

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']


def setup(app):
    app.add_css_file('theme_overrides.css')
