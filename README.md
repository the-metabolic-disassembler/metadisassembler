# The Metabolic Disassembler

The Metabolic Disassembler is a Python package to automatically predict a combination of biosynthetic units in a natural product. This software would help to reveal the starting materials of the target natural product.  


## Basic Usage
A usage example by Jupyter Notebook can be seen [here](https://github.com/the-metabolic-disassembler/metadisassembler/blob/master/jupyter_usecase/basic_usage.ipynb).   

<a href="https://github.com/the-metabolic-disassembler/metadisassembler/blob/master/jupyter_usecase/basic_usage.ipynb"><img src="https://raw.githubusercontent.com/the-metabolic-disassembler/metadisassembler/master/jupyter_usecase/jupyter01.png" width="400"></a>


## Command Line Usage

```bash
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
```

## Requirements

- Python (3.6)
  - [RDKit](https://www.rdkit.org) (version 2018.09.1.0)
  - [NetworkX](https://networkx.github.io/documentation/stable/) (version 2.2)
  - [CairoSVG](https://cairosvg.org) (version 2.2.1)
  - [Pillow (PIL)](https://pillow.readthedocs.io/en/stable/)
  - [Pandas](https://pandas.pydata.org)
  - [Matplotlib](https://matplotlib.org)

## Installation

Install The Metabolic Disassembler with `pip`.  

```bash
$ pip install metadisassembler
```

## Google Colab

You can try the process from installation to the basic usage in Google Colab. [![colab-logo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/the-metabolic-disassembler/metadisassembler/blob/master/jupyter_usecase/basic_usage_in_colab.ipynb)

<a href="https://colab.research.google.com/github/the-metabolic-disassembler/metadisassembler/blob/master/jupyter_usecase/basic_usage_in_colab.ipynb"><img src="https://raw.githubusercontent.com/the-metabolic-disassembler/metadisassembler/master/jupyter_usecase/colab01.png" width="400"></a>
