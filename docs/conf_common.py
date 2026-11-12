"""Shared Sphinx settings for independently maintained English and Chinese docs."""

from pathlib import Path
import sys

# -- Path setup --------------------------------------------------------------

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# -- Project information -----------------------------------------------------

project = "NEP-kappa"
copyright = "2026, Shixian Liu, Fei Yin"
author = "Shixian Liu, Fei Yin"
release = "2.0.1"

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.graphviz",
]

templates_path = ["_templates"]
exclude_patterns = []

source_suffix = ".rst"
master_doc = "index"

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_logo = "_static/NEP-kappa-logo-1200x750.png"

html_theme_options = {
    "navigation_depth": 3,
    "collapse_navigation": True,
    "sticky_navigation": True,
    "titles_only": False,
    "logo_only": True,
}

html_short_title = "NEP-kappa"

# Vector diagrams remain sharp when zoomed and preserve Unicode labels.
graphviz_output_format = "svg"

# -- Autodoc options ---------------------------------------------------------

autodoc_member_order = "bysource"
napoleon_google_docstring = True
napoleon_numpy_docstring = True
