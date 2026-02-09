# docs/source/conf.py
from __future__ import annotations

import os
import sys

# -- Path setup --------------------------------------------------------------
# If your package uses a src/ layout, add it so autodoc can import flowsig.
sys.path.insert(0, os.path.abspath("../../src"))

# -- Project information -----------------------------------------------------
project = "flowsig"
author = "Axel A. Almet"
copyright = "2024, Axel A. Almet"

# Pull version from installed package metadata if available.
# This works well on RTD when you `pip install -e . --no-deps`.
try:
    from importlib.metadata import version as _version

    release = _version("flowsig")
except Exception:
    release = "unknown"

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",            # Google/Numpy style docstrings
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",       # show type hints nicely
    "myst_parser",                    # allow Markdown pages
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Autodoc / autosummary
autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}
autodoc_typehints = "description"

# IMPORTANT: your stack is heavy. Mock imports so docs build without installing
# TensorFlow/scanpy/squidpy/geopandas/etc.
autodoc_mock_imports = [
    "tensorflow",
    "tensorflow_probability",
    "tf_keras",
    "scanpy",
    "squidpy",
    "anndata",
    "geopandas",
    "pyliger",
    "cnmf",
    "dm_tree",
    "spatial_factorization",
    "causaldag",
    "graphical_models",
    "graphical_model_learning",
]

# MyST (Markdown) settings (safe defaults)
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
]
myst_heading_anchors = 3

# Intersphinx links (optional but useful)
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy", None),
}

# -- Options for HTML output -------------------------------------------------
html_theme = "shibuya"
