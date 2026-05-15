[![PyPI](https://img.shields.io/pypi/v/physquirrel)](https://pypi.org/project/physquirrel/)
[![License](https://img.shields.io/github/license/nholtgrefe/squirrel)](https://github.com/nholtgrefe/squirrel/blob/main/LICENSE)
[![Docs](https://img.shields.io/badge/docs-dev-blue)](https://github.com/nholtgrefe/squirrel/blob/main/DOCS.md)
[![MBE DOI](https://img.shields.io/badge/MBE-10.1093%2Fmolbev%2Fmsaf067-blue)](https://doi.org/10.1093/molbev/msaf067)

# physquirrel

physquirrel is the Python library for **Squirrel**—an algorithm for reconstructing semi-directed
phylogenetic level-1 networks from quarnets and/or sequence alignments. It uses
[networkx](https://networkx.org/) for network representations and [numpy](https://numpy.org/),
[scipy](https://scipy.org/), and [matplotlib](https://matplotlib.org/) for computation and
visualization.

<br>

## Key Features

- **δ-heuristic**: construct quarnets (4-leaf subnetworks) from multiple sequence alignments in `.fasta` or `.nexus` format
- **Squirrel algorithm**: reconstruct semi-directed phylogenetic level-1 networks from quarnets
- **Visualization**: basic plotting of phylogenetic networks
- **eNewick export**: serialize phylogenetic trees and networks in `eNewick` format

## Installation

```bash
pip install physquirrel
```

Runtime dependencies (`numpy`, `networkx`, `matplotlib`, `scipy`) are installed automatically.

## Documentation

For the full manual and API reference, visit the **[physquirrel docs](https://github.com/nholtgrefe/squirrel/blob/main/DOCS.md)**.

## Citation

If you use physquirrel, please cite:

> Niels Holtgrefe, Katharina T. Huber, Leo van Iersel, Mark Jones, Samuel Martin, and Vincent Moulton.
> **Squirrel: Reconstructing semi-directed phylogenetic level-1 networks from four-leaved networks or sequence alignments.**
> *Molecular Biology and Evolution*, 42(4):msaf067, 2025. doi: [10.1093/molbev/msaf067](https://doi.org/10.1093/molbev/msaf067)

## See also

For the graphical user interface developed for the paper, please go to [`gui/`](https://github.com/nholtgrefe/squirrel/tree/main/gui).

For the experimental materials corresponding to the paper, please go to [`experiments/`](https://github.com/nholtgreve/squirrel/tree/main/experiments).
