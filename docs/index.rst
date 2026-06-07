

.. image:: _static/squirrel_icon.svg
   :alt: physquirrel logo
   :align: right
   :width: 150px


**Welcome to the physquirrel docs!**

physquirrel is the Python implementation of the **Squirrel** algorithm for
reconstructing level-1 phylogenetic networks from quartet data and/or sequence alignments. It works on
networks, MSAs, and distance matrices represented by `PhyloZoo
<https://github.com/nholtgrefe/phylozoo>`_, infers 4-leaf subnetworks via the
delta heuristic, and then reconstructs semi-directed or rooted networks.

Documentation Overview
-----------------------

New to physquirrel? Start with the :doc:`Quickstart <quickstart>`.
For a complete walkthrough of all input types and options, see the
:doc:`Manual <manual>`. For detailed installation instructions, see the :doc:`Installation Guide <installation>`.
See the :doc:`API Reference <api/index>` for detailed documentation of all public
symbols, organized by module. 

.. toctree::
   :maxdepth: 2
   :caption: Contents

   quickstart
   installation
   manual
   api/index
   changelog

Indices and tables
------------------

* :ref:`genindex` — index of all functions and classes.
* :ref:`modindex` — index of all modules.
* :ref:`search` — search the documentation.

Citation
--------

If you use this package in research, please cite:

   Niels Holtgrefe, Katharina T. Huber, Leo van Iersel, Mark Jones, Samuel Martin, and Vincent Moulton.
   **Squirrel: Reconstructing semi-directed phylogenetic level-1 networks from four-leaved networks or sequence alignments.**
   *Molecular Biology and Evolution*, 42(4):msaf067, 2025.
   doi: `10.1093/molbev/msaf067 <https://doi.org/10.1093/molbev/msaf067>`_



.. image:: _static/squirrel_cover.svg
   :alt: physquirrel logo
   :align: left
