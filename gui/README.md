# Graphical User-Interface (GUI) for Squirrel

The **Squirrel-GUI** is a Graphical User-Interface for phylogenetic network analysis, employing the Squirrel algorithm to reconstruct semi-directed phylogenetic level-1 networks from quarnets or sequence alignments.

## Installation
The GUI is available for Windows and Linux [here](https://drive.google.com/drive/folders/16ZQ2qM0v4-DgKY3A1av1ms9NeEKSDe1C?usp=sharing). To install and run, simply download the `.exe` file for your operating system and open/execute it once downloaded. (If preffered, experienced users can also run the GUI from the [source files](https://github.com/nholtgrefe/squirrel/tree/main/gui/src) or use our Python package [`physquirrel`](https://github.com/nholtgrefe/squirrel/tree/main/physquirrel).)

## Features

### File Window
- `Select file`: Allows the user to select a file out of one of the two following options:
	1.  `.fasta` or `.nexus` file containing a multiple sequence alignment;
	2.  `.txt` file containing a dense set of tf-quarnets, where each line contains one tf-quarnet in the folowing format:
		- `SQ: a b c d 1.0` for a quarnet on leaves $\{a,b,c,d\}$ with a split $ab|cd$ and weight 1.0;
		- `4C: a b c d 1.0` for a quarnet on leaves $\{a,b,c,d\}$ with a four-cycle $a,b,c,d$, the leaf $a$ below the reticulation and weight 1.0.
- `File information`: Displays information on the chosen file, such as a list of taxa and the length of the sequence alignment.

<img src="https://github.com/user-attachments/assets/2ab886b6-00c6-482d-9fb4-7d5ccd02b6f1" alt="Sample Image" width="700" >

### Algorithm Window
- `Reconstruct networks`: Runs the **Squirrel** algorithm (and possibly the **$\delta$-heuristic**) to reconstruct a phylogenetic level-1 network from the set of tf-quarnets or the sequence alignment.
- `Algorithm settings`: Allows the user to:
	- specify an outgroup (which results in rooted networks)
	- change the $\lambda$-value for the $\delta$-heuristic,
	- change the maximum number of leaves for which the Travelling Salesman Problem is solved optimally
	- choose whether the quarnet-weights should be used.
- `Algorithm status`: Window to show the progress while the algorithm is running.
![gui02](https://github.com/user-attachments/assets/f3f4480e-167c-4df4-9abe-f7147e8d76a7)

### Network Window
- `Reconstructed networks`: Table with the reconstructed networks and corresponding (weighted) consistency scores. The table will continuously be updated while the algorithm is running. The table also allows the user to sort the networks according to their scores.
- `eNewick save options`: Allows the user to either save the selected network or all networks into a `.txt` file using the `eNewick` format. In case no outgroup was specified, the semi-directed networks are rooted randomly to obtain an `eNewick` string.
- `Network visualizer`: Basic visualization window to depict the selected network and its corresponding `eNewick` string. The visualization is not optimized for larger networks, in which case we recommend an external program (such as [Dendroscope 3](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/dendroscope/)) to view the networks in more detail.
![gui03](https://github.com/user-attachments/assets/df1775ad-b442-448f-be77-5665ee5634ca)

## Citation
If you use the **Squirrel** GUI, please cite the corresponding paper:

*Squirrel: Reconstructing semi-directed phylogenetic level-1 networks from four-leaved networks and sequence alignments* by Niels Holtgrefe, Katharina T. Huber, Leo van Iersel, Mark Jones, Samuel Martin, and Vincent Moulton.
