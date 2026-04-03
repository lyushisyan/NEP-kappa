from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

project = "NEP-kappa"
author = "Shixian Liu, Fei Yin"
release = "0.1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
]

templates_path = ["_templates"]
exclude_patterns = []

language = "en"

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

autodoc_member_order = "bysource"
autodoc_mock_imports = [
    "ase",
    "calorine",
    "hiphive",
    "trainstation",
    "phonopy",
    "phono3py",
    "h5py",
    "matplotlib",
    "seekpath",
]
