Installation
============

Installing physquirrel
-----------------------

physquirrel is a Python package that runs on `Python <https://www.python.org/>`_
(>= 3.10). Choose one of:

* **From PyPI** — stable release:

  .. code-block:: bash

     pip install physquirrel

* **From source** — editable install, recommended for development:

  .. code-block:: bash

     git clone https://github.com/nholtgrefe/squirrel
     cd squirrel
     pip install -e .

Optional dependency groups
~~~~~~~~~~~~~~~~~~~~~~~~~~

Visualization extras (via PhyloZoo's plotting utilities):

.. code-block:: bash

  pip install physquirrel[viz]

Graphviz extras (via PhyloZoo's plotting utilities, for Graphviz-based rendering):

.. code-block:: bash

  pip install physquirrel[graphviz]

Development and testing tools:

.. code-block:: bash

   pip install physquirrel[dev]

Documentation dependencies:

.. code-block:: bash

   pip install physquirrel[docs]

Requirements
^^^^^^^^^^^^

The following are required and installed automatically with pip:

* `phylozoo <https://github.com/nholtgrefe/phylozoo>`_ >= 0.1.2 — core datatypes
  for phylogenetic networks, splits, quartets, distance matrices, and MSAs.
* `numpy <https://numpy.org/>`_ >= 1.20 — numerical arrays and distance matrix
  operations.
* `networkx <https://networkx.org/>`_ >= 3.0 — graph algorithms used in quartet
  joining and cycle resolution.
* `numba <https://numba.pydata.org/>`_ >= 0.56 — JIT-accelerated Held-Karp TSP
  solver used in cycle resolution.

Verifying Installation
-----------------------

To verify that physquirrel is installed correctly, import it and print the version.
The latest version is |version|.

.. code-block:: python

   >>> import physquirrel as sq
   >>> print(sq.__version__)
   x.y.z  # your installed version

Building Documentation
-----------------------

To build the documentation locally, install the optional documentation
dependencies and run sphinx-build from the repository root:

.. code-block:: bash

   pip install -e ".[docs]"
   sphinx-build -b html docs docs/_build/html

Open ``docs/_build/html/index.html`` in a browser.

Troubleshooting
---------------

**Slow first call**: The Held-Karp TSP kernel (used inside
:func:`~physquirrel.algorithms.tsp.optimal_tsp_tour`) is JIT-compiled by Numba
on first invocation. This is expected — subsequent calls in the same Python
session are fast.
