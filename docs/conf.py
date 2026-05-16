"""Sphinx configuration for physquirrel documentation."""

from __future__ import annotations

import importlib.metadata

try:
    release = importlib.metadata.version("physquirrel")
except importlib.metadata.PackageNotFoundError:
    release = "2.0.0"

project = "physquirrel"
author = "Niels Holtgrefe"
copyright = "2024–2026, Niels Holtgrefe"
version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
]

autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}

templates_path = ["_templates"]
exclude_patterns: list[str] = ["_build"]

html_theme = "pydata_sphinx_theme"
html_title = "physquirrel"
html_theme_options = {
    "logo": {
        "image_light": "_static/squirrel_icon.svg",
        "image_dark": "_static/squirrel_icon.svg",
        "text": "physquirrel",
        "alt_text": "physquirrel",
    },
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/nholtgrefe/squirrel",
            "icon": "fa-brands fa-github",
            "type": "fontawesome",
        },
    ],
}
html_static_path = ["_static"]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "networkx": ("https://networkx.org/documentation/stable/", None),
    "phylozoo": ("https://nholtgrefe.github.io/phylozoo/", None),
}

nitpick_ignore = [
    ("py:class", "phylozoo.core.network.sdnetwork.SemiDirectedPhyNetwork"),
    ("py:class", "phylozoo.core.network.dnetwork.DirectedPhyNetwork"),
    ("py:class", "phylozoo.core.quartet.qprofile.QuartetProfile"),
    ("py:class", "phylozoo.core.quartet.qprofileset.QuartetProfileSet"),
    ("py:class", "phylozoo.core.primitives.circular_ordering.CircularOrdering"),
    ("py:class", "phylozoo.core.primitives.circular_ordering.CircularSetOrdering"),
    ("py:class", "phylozoo.core.primitives.partition.Partition"),
    ("py:class", "phylozoo.utils.parallel.ParallelConfig"),
]

suppress_warnings = ["ref.python"]
