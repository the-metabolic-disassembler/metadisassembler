The Metabolic Disassembler
==========================

The Metabolic Disassembler is a Python package to automatically predict
a combination of biosynthetic building blocks in a metabolic compound.
This software would help to reveal the basic metabolites constructing
the target product.

Installation
------------

Install The Metabolic Disassembler with ``pip``.

.. code:: bash

   $ pip install metadisassembler

Requirements
------------

-  Python (3.6)

   -  `RDKit <https://www.rdkit.org>`__
   -  `NetworkX <https://networkx.github.io/documentation/stable/>`__
   -  `CairoSVG <https://cairosvg.org>`__
   -  `Pillow (PIL) <https://pillow.readthedocs.io/en/stable/>`__
   -  `Pandas <https://pandas.pydata.org>`__
   -  `Matplotlib <https://matplotlib.org>`__

Command Line Usage
------------------

.. code:: bash

   % metadisassembler -h
   usage: metadisassembler [-h] [-t TIME] [--hide] [-c] [-p] query

   positional arguments:
     query                 MDL_Molfile, SMILES, InChI, KEGG_COMPOUND_ID,
                           KNApSAcK_COMPOUND_ID

   optional arguments:
     -h, --help            show this help message and exit
     -t TIME, --time TIME  set a time limit [s] [default: 300]
     --hide                hide stereochemistry [default: False]
     -c, --color           output color allocation information [default: False]
     -p, --pro             you never miss any combination (NOT RECOMMENDED).
                           [default: False]

Basic Usage
-----------

An example notebook is available
`here <https://github.com/the-metabolic-disassembler/metadisassembler/blob/master/jupyter_usecase/basic_usage.ipynb>`__.

You can try it in Google Colab. |colab-logo|

.. |colab-logo| image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/github/the-metabolic-disassembler/metadisassembler/blob/master/jupyter_usecase/basic_usage_in_colab.ipynb
