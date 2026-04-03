# Configuration file for the Sphinx documentation builder.

import os
import sys

# -- Path setup --------------------------------------------------------------

sys.path.insert(0, os.path.abspath("../.."))

# -- Project information -----------------------------------------------------

project = "NEP-kappa"
copyright = "2026, Shixian Liu, Fei Yin"
author = "Shixian Liu, Fei Yin"
release = "0.1.0"

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

templates_path = ["_templates"]
exclude_patterns = []

source_suffix = ".rst"
master_doc = "index"
language = "en"

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

html_theme_options = {
    "navigation_depth": 4,
    "collapse_navigation": False,
    "sticky_navigation": True,
    "titles_only": False,
}

html_title = "NEP-kappa Documentation"
html_short_title = "NEP-kappa"

# -- Autodoc options ---------------------------------------------------------

autodoc_member_order = "bysource"
napoleon_google_docstring = True
napoleon_numpy_docstring = True