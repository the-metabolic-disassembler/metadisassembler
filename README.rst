The Metabolic Disassembler
==========================

The Metabolic Disassembler is a Python package to automatically predict
a combination of biosynthetic units in a natural product. This software
would help to reveal the starting materials of the target natural
product.

Basic Usage
-----------

A usage example by Jupyter Notebook can be seen
`here <https://github.com/the-metabolic-disassembler/metadisassembler/blob/master/jupyter_usecase/basic_usage.ipynb>`__.

Command Line Usage
------------------

.. code:: bash

   % metadisassembler -h
   usage: metadisassembler [-h] [-t TIME] [--hide] [-c] query

   positional arguments:
     query                 MDL_Molfile, SMILES, InChI, KEGG_COMPOUND_ID,
                           KNApSAcK_COMPOUND_ID

   optional arguments:
     -h, --help            show this help message and exit
     -t TIME, --time TIME  set a time limit [s] [default: 300]
     --hide                hide stereochemistry [default: False]
     -c, --color           output color allocation information [default: False]

Requirements
------------

-  Python (3.6)

   -  `RDKit <https://www.rdkit.org>`__ (version 2018.09.1.0)
   -  `NetworkX <https://networkx.github.io/documentation/stable/>`__
      (version 2.2)
   -  `CairoSVG <https://cairosvg.org>`__ (version 2.2.1)
   -  `Pillow (PIL) <https://pillow.readthedocs.io/en/stable/>`__
   -  `Pandas <https://pandas.pydata.org>`__
   -  `Matplotlib <https://matplotlib.org>`__

Installation
------------

Install The Metabolic Disassembler with ``pip``.

.. code:: bash

   $ pip install metadisassembler

Google Colab
------------

You can try the process from installation to the basic usage in Google
Colab. |colab-logo|

.. |colab-logo| image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/github/the-metabolic-disassembler/metadisassembler/blob/master/jupyter_usecase/basic_usage_in_colab.ipynb
