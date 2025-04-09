# `physquirrel`
`physquirrel` is a Python package for phylogenetic network analysis, with a focus on reconstructing semi-directed phylogenetic level-1 networks from quarnets and/or sequence alignments. The package provides tools to build and visualize phylogenetic networks, leveraging the **Squirrel** algorithm for efficient network reconstruction.


## List of important features
- **$\delta$-heuristic** to construct quarnets (4-leaf subnetworks) from a multiple sequence alignment (in `.fasta` or `.nexus` format)
- **Squirrel** algorithm to construct semi-directed phylogenetic level-1 networks from quarnets
- Basic visualization of networks
- Exporting phylogenetic trees and networks in `eNewick` format
- Methods to extract information from a network (e.g. its set of splits, its displayed quarnets)

## Installation
If you have an up-to-date version of [Python](https://www.python.org/downloads/) installed on your device, the standard package manager `pip` should come pre-installed. Then, you can install `physquirrel` from [PyPI](https://pypi.org/project/physquirrel/) by simply using the following command in a terminal:

```
python -m pip install physquirrel
```

  
## Example usage

### Importing the package
To get started with `physquirrel`, open a Python shell and import the package with:

```
import physquirrel as psq
```

### Creating a set of quarnets
Use the $\delta$-heuristic to create a dense set of tf-quarnets from a multiple sequence alignment as follows. (To test out the package, we recommend using the small sequence alignment in the file [`hiv.fasta`](https://github.com/nholtgrefe/squirrel/blob/main/data/hiv/hiv.fasta).
:
```
msa = psq.MSA('path/to/msa/file.fasta')
Q = msa.delta_heuristic()
```

Alternatively, the dense set of tf-quarnets can also be loaded directly from a `.txt` file as follows:
```
Q = psq.DenseQuarnetSet('path/to/quarnet/file.txt')
```
This method assumes that the `.txt` file contains one line per tf-quarnet. The quarnets need to be one of the following two types:
1. `SQ: a b c d` for a quarnet on leaves $\{a,b,c,d\}$ with a split $ab|cd$.
2. `4C: a b c d` for a quarnet on leaves $\{a,b,c,d\}$ with a four-cycle $a,b,c,d$ and the leaf $a$ below the reticulation.

To give the quarnets a weight, simply add the weight to the end of the string (e.g. `SQ: a b c d 0.5`).

### Reconstructing a network
To create a phylogenetic network from the dense set of tf-quarnets, run the Squirrel algorithm. (If some outgroup `taxon_name` is known, this can be specified by passing `outgroup='taxon_name'` to the method.)
```
N = Q.squirrel()
```
To view the resulting network and show its `eNewick` string (with an arbitrary rooting if no outgroup was specified), run:
```
N.visualize()
N.create_enewick()
```
        
For a complete overview of different methods and extra parameter options, please check the method descriptions in the [source code](https://github.com/nholtgrefe/squirrel/tree/main/physquirrel/src/physquirrel) of `physquirrel`.


## Citation
If you use `physquirrel`, please cite the corresponding paper:

> **Squirrel: Reconstructing semi-directed phylogenetic level-1 networks from four-leaved networks or sequence alignments**.
> *Niels Holtgrefe, Katharina T. Huber, Leo van Iersel, Mark Jones, Samuel Martin, and Vincent Moulton.*
> Molecular Biology and Evolution, 42(4):msaf067, 2025. doi: [10.1093/molbev/msaf067](https://doi.org/10.1093/molbev/msaf067)
