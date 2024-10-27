# Graphical User-Interface (GUI) for Squirrel

This page describes the Graphical User-Interface (GUI) that comes with the **Squirrel** algorithm.

## Installation
The GUI is available for Windows and Linux [here](https://drive.google.com/drive/folders/16ZQ2qM0v4-DgKY3A1av1ms9NeEKSDe1C?usp=sharing). To install and run, simply download the `.exe` file for your operating system and open it once downloaded.

## Features

### File Window
- `Select file`: Allows the user to select a file out of one of the two following options:
	1.  `.fasta` or `.nexus` file containing a multiple sequence alignment;
	2.  `.txt` file containing a dense set of tf-quarnets, where each line contains one tf-quarnet in the folowing format:
		- `SQ: a b c d 1.0` for a quarnet on leaves $\{a,b,c,d\}$ with a split $ab|cd$ and weight 1.0;
		- `4C: a b c d 1.0` for a quarnet on leaves $\{a,b,c,d\}$ with a four-cycle $a,b,c,d$, the leaf $a$ below the reticulation and weight 1.0.
- `File information`: Displays information on the chosen file, such as a list of taxa and the length of the sequence alignment.

![gui1](https://github.com/user-attachments/assets/6daa958e-a460-4260-a786-4a01302a30c1)

### Algorithm Window
- `Reconstruct networks`: Runs the **Squirrel** algorithm (and possibly the **$\delta$-heuristic**) to reconstruct a phylogenetic level-1 network from the set of tf-quarnets or the sequence alignment.
- `Algorithm settings`: Allows the user to:
	- specify an outgroup (which results in rooted networks)
	- change the $\lambda$-value for the $\delta$-heuristic,
	- change the maximum number of leaves for which the Travelling Salesman Problem is solved optimally
	- choose whether the quarnet-weights should be used.
- `Algorithm status`: Window to show the progress while the algorithm is running.
![gui2](https://github.com/user-attachments/assets/26d2847b-608e-4d9d-becf-412baf864280)

### Network Window
- `Reconstructed networks`: Table with the reconstructed networks and corresponding (weighted) consistency scores. The table will continuously be updated while the algorithm is running. The table also allows the user to sort the networks according to their scores.
- `eNewick save options`: Allows the user to either save the selected network or all networks into a `.txt` file using the `eNewick` format. In case no outgroup was specified, the semi-directed networks are rooted randomly to obtain an `eNewick` string.
- `Network visualizer`: Basic visualization window to depict the selected network and its corresponding `eNewick` string. The visualization is not optimized for larger networks, in which case we recommend an external program (such as [Dendroscope 3](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/dendroscope/)) to view the networks in more detail.
![gui3](https://github.com/user-attachments/assets/ee492e3d-e3df-48e9-8fd6-0b5bd1236d9e)

## Citation
If you use the **Squirrel** GUI, please cite the corresponding paper:

*Squirrel: Reconstructing semi-directed phylogenetic level-1 networks from four-leaved networks and sequence alignments* by Niels Holtgrefe, Katharina T. Huber, Leo van Iersel, Mark Jones, Samuel Martin, and Vincent Moulton.
