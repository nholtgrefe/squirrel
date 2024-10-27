# Squirrel
This is the repository corresponding to the paper: 'Squirrel: Reconstructing semi-directed phylogenetic level-1 networks from four-leaved networks and sequence alignments' by Niels Holtgrefe, Katharina T. Huber, Leo van Iersel, Mark Jones, Samuel Martin, and Vincent Moulton.

## Running the Squirrel algorithm
The **Squirrel** algorithm can either be run with a graphical user-interface (GUI) for inexperienced programmers, or as part of our Python package `physquirrel` which allows for more functionality (and is therefore recommended if the user is experienced with Python). For an installation guide and a more detailed overview of the features, please visit the `README.md` files in the respective folders:
- [gui](https://github.com/nholtgrefe/squirrel/tree/main/gui), for the graphical user-interface;
- [physquirrel], for the corresponding Python package.

To test out the GUI or the Python package, we recommend using the small sequence alignment in the `hiv.fasta` [file].

## Supplementary files
The biological datasets that appeared in the paper are in the folder [data], together with the networks and results of the simulation experiments. The corresponding Python scripts are in the folder [scripts]. For a more detailed explanation of these files, see the respective `README.md` files in the folders.

## Citation
If you use any of the code in this repository, please cite the corresponding paper: 'Squirrel: Reconstructing semi-directed phylogenetic level-1 networks from four-leaved networks and sequence alignments' by Niels Holtgrefe, Katharina T. Huber, Leo van Iersel, Mark Jones, Samuel Martin, and Vincent Moulton.

Or in `bibtex` format:
```
@misc{holtgrefe2024squirrel,
	 title={Squirrel: Reconstructing semi-directed phylogenetic level-1 networks from four-leaved networks and sequence alignments}, 
	 author={Holtgrefe, Niels and Huber, Katharina T, and van Iersel, Leo and Jones, Mark and Martin, Samuel and Moulton, Vincent},
	 year={2024}
	 }
```
