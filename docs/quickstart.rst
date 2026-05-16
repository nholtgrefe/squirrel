Quickstart
==========

This page shows the minimal steps to go from a sequence alignment to a
phylogenetic network in physquirrel. For more details, 
see the :doc:`Manual <manual>` and the :doc:`Installation <installation>` page.

Install
-------

.. code-block:: bash

   pip install physquirrel

Run Squirrel
------------

.. code-block:: python

   import physquirrel as psq

   # Load a multiple sequence alignment (.fasta, .fa, .fas, .nexus, .nex)
   msa = psq.MSA.load("alignment.fasta")

   # Reconstruct a semi-directed level-1 network
   network = psq.squirrel_from_msa(msa)

   # Print the eNewick string
   print(network.to_string())

   # Save to file
   network.save("network.nwk")

That is all that is needed. By default the network is unrooted
(:class:`~phylozoo.SemiDirectedPhyNetwork`). To root it, pass an outgroup taxon
name:

.. code-block:: python

   rooted = psq.squirrel_from_msa(msa, outgroup="Taxon1")
   rooted.save("network_rooted.nwk")

Visualize (optional)
--------------------

Plotting relies on PhyloZoo and requires the ``viz`` extra dependencies (``pip install physquirrel[viz]`` or
``pip install physquirrel[graphviz]``); see the :doc:`Installation Guide <installation>`.

.. code-block:: python

   from phylozoo.viz import plot

   plot(network)

Next steps
----------

* :doc:`manual` — full workflow, all input types, and I/O options
* :doc:`installation` — optional dependencies and building from source
* :doc:`api/index` — complete API reference
