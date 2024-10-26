# `physquirrel`
`physquirrel` is a Python package for phylogenetic network analysis, with a focus on reconstructing semi-directed phylogenetic level-1 networks from quarnets and/or sequence alignments. The package provides tools to build and visualize phylogenetic networks, leveraging the **Squirrel** algorithm for efficient network reconstruction.


## List of important features
- **$\delta$-heuristic** to construct quarnets (4-leaf subnetworks) from a multiple sequence alignment (in `.fasta` or `.nexus` format)
- **Squirrel** algorithm to construct semi-directed phylogenetic level-1 networks from quarnets
- Visualization of networks
- Exporting phylogenetic trees and networks in `eNewick` format
- Methods to extract information from a network (e.g. its set of splits, its displayed quarnets)

## Installation
If you have an up-to-date version of [Python](https://www.python.org/downloads/) installed on your device, the standard package manager `pip` should come pre-installed. Then, you can install `physquirrel` from [PyPI](https://pypi.org/) by simply using the following command in a terminal:

`python -m pip install physquirrel`

  
## Example usage
To get started with `physquirrel`, open a Python shell and import the package with:

```
import physquirrel as psq
```
Then, use the $\delta$-heuristic to create a dense set of tf-quarnets from a multiple sequence alignment as follows:
```
msa = psq.MSA('path/to/msa/file.fasta')
Q = msa.delta_heuristic()
```
To create a network from the dense set of tf-quarnets, run the Squirrel algorithm note that an outgroup can optionally be provided with the variable `outgroup`):
```
N = Q.squirrel()
```
To view the network and print its `eNewick` string (with an arbitrary rooting), run:
```
N.visualize()
eNewick = N.create_enewick()
print(eNewick)
```

Alternatively, the dense set of tf-quarnets can also be loaded directly from a `.txt` file as follows:
```
Q = psq.DenseQuarnetSet('path/to/quarnet/file.txt')
```
This method assumes that the `.txt` file contains one line per tf-quarnet. The quarnets need to be one of the following two types:
1. `SQ: a b c d` for a quarnet on leaves $\{a,b,c,d\}$ with a split $ab|cd$.
2. `4C: a b c d` for a quarnet on leaves $\{a,b,c,d\}$ with a four-cycle $a,b,c,d$ and the leaf $a$ below the reticulation.
        
For a complete overview of different methods and extra parameter options, please check the method descriptions in the [source code](https://github.com/nholtgrefe/squirrel/tree/main/physquirrel/src/physquirrel) of `physquirrel`.


## Citation
If you use `physquirrel`, please cite the corresponding paper:

*Squirrel: Reconstructing semi-directed phylogenetic level-1 networks from four-leaved networks and sequence alignments* by Niels Holtgrefe, Katharina T. Huber, Leo van Iersel, Mark Jones, Samuel Martin, and Vincent Moulton.
