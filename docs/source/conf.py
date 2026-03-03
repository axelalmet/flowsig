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
    "sphinx.ext.napoleon",            # Google/Numpy style docstrings
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",       # show type hints nicely
    "myst_parser",                    # allow Markdown pages
    "autoapi.extension"
]     

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- AutoAPI configuration ---------------------------------------------------
autoapi_type = "python"
autoapi_dirs = ["../../src/flowsig"]   # path from docs/source to your package
autoapi_root = "autoapi"
autoapi_packages = ["flowsig"]
autoapi_options = [
    "members",
    "undoc-members",
    "show-inheritance",
    "show-module-summary",
    "imported-members",
]
autoapi_ignore = [
    "*/_version*",
    "*/tests/*",
    "*example_script*",          # ignore example scripts
    "*-checkpoint*",             # ignore Jupyter checkpoint artifacts
]
autoapi_keep_files = True              # keep generated .rst for inspection
autoapi_python_class_content = "both"  # show class docstring + __init__

# Autodoc / autosummary
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
    "mofaflex",
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

# -- AutoAPI skip handler (replaces the old autodoc one) ---------------------
def _autoapi_skip_member(app, what, name, obj, skip, options):
    """Keep public API visible, hide private implementation details."""
    short_name = name.split(".")[-1]
    # Skip private functions/classes/attributes (but NOT private modules)
    if what in ("function", "method", "attribute", "property", "class"):
        if short_name.startswith("_") and not short_name.startswith("__"):
            return True
    # Skip dunder methods except __init__
    if short_name.startswith("__") and short_name not in ("__init__",):
        return True
    return skip

def setup(app):
    app.connect("autoapi-skip-member", _autoapi_skip_member)