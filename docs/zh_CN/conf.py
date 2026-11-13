"""简体中文文档的独立构建配置，共用英文版的主题和资源。"""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from conf_common import *

language = "zh_CN"
html_title = "NEP-kappa 中文文档"
html_static_path = ["../source/_static"]
html_logo = "../source/_static/NEP-kappa-logo-1200x750.png"
html_css_files = ["chinese.css"]
