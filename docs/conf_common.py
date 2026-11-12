"""Shared Sphinx settings for independently maintained English and Chinese docs."""

from pathlib import Path
import os
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

templates_path = [str(Path(__file__).resolve().parent / "_templates")]
exclude_patterns = []

source_suffix = ".rst"
master_doc = "index"

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_css_files = ["language-switch.css"]
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


def language_links(app, pagename, templatename, context, doctree):
    """Link matching pages across local builds and Read the Docs translations."""
    current = "zh-cn" if app.config.language == "zh_CN" else "en"
    target = app.builder.get_target_uri(pagename)
    hosted = os.environ.get("READTHEDOCS") == "True"
    version = os.environ.get("READTHEDOCS_VERSION", "latest")
    links = []
    for code, label, directory in (("en", "English", "en"), ("zh-cn", "中文", "zh_CN")):
        if hosted:
            url = f"https://nep-kappa.readthedocs.io/{code}/{version}/{target}"
        else:
            destination = Path(app.outdir).parent / directory / target
            url = Path(os.path.relpath(destination, (Path(app.outdir) / target).parent)).as_posix()
        links.append({"label": label, "code": code, "url": url, "current": code == current})
    context["language_links"] = links
    context["language_label"] = "语言" if current == "zh-cn" else "Language"


def setup(app):
    app.connect("html-page-context", language_links)
