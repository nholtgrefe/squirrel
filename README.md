# Squirrel
This is the repository corresponding to the paper: 
> **Squirrel: Reconstructing semi-directed phylogenetic level-1 networks from four-leaved networks or sequence alignments**.
> *Niels Holtgrefe, Katharina T. Huber, Leo van Iersel, Mark Jones, Samuel Martin, and Vincent Moulton.* 2024. bioRxiv: [10.1101/2024.11.01.621567](https://doi.org/10.1101/2024.11.01.621567)

## Using the Squirrel algorithm
The **Squirrel** algorithm is easy to install and can either be run with a graphical user-interface (GUI) for inexperienced programmers, or as part of our Python package `physquirrel` which allows for more functionality (and is therefore recommended if the user is experienced with Python). For an installation guide and a more detailed overview of the features (including screenshots for the GUI), please click on the links below to visit the `README.md` files in the respective folders:
- [gui](https://github.com/nholtgrefe/squirrel/tree/main/gui), for the graphical user-interface (Windows and Linux only);
- [physquirrel](https://github.com/nholtgrefe/squirrel/tree/main/physquirrel), for the corresponding Python package.

To test out the GUI or the Python package, we recommend using the small sequence alignment in the file [`hiv.fasta`](https://github.com/nholtgrefe/squirrel/blob/main/data/hiv/hiv.fasta).
## Supplementary data and scripts
The biological data sets, simulation results and used Python scripts from the paper are all in the folder [data](https://github.com/nholtgrefe/squirrel/tree/main/data). For a more detailed explanation of these files (including sources for the biological data), see the `README.md` file in the folder.

## Citation
If you use any of the code in this repository, please cite the corresponding paper: 'Squirrel: Reconstructing semi-directed phylogenetic level-1 networks from four-leaved networks and sequence alignments' by Niels Holtgrefe, Katharina T. Huber, Leo van Iersel, Mark Jones, Samuel Martin, and Vincent Moulton. 2024. bioRxiv: [10.1101/2024.11.01.621567](https://doi.org/10.1101/2024.11.01.621567)

Or in `bibtex` format:
```
@misc{holtgrefe2024squirrel,
	 title={Squirrel: Reconstructing semi-directed phylogenetic level-1 networks from four-leaved networks or sequence alignments}, 
	 author={Holtgrefe, Niels and Huber, Katharina T, and van Iersel, Leo and Jones, Mark and Martin, Samuel and Moulton, Vincent},
	 year={2024},
         archivePrefix={bioRxiv},
	 eprint={10.1101/2024.11.01.621567}
	 }
```
