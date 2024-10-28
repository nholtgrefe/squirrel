# Supplementary data and scripts
This folder contains the biological datasets, the simulated networks and results of the simulation experiments from the paper:
'Squirrel: Reconstructing semi-directed phylogenetic level-1 networks from four-leaved networks and sequence alignments' by Niels Holtgrefe, Katharina T. Huber, Leo van Iersel, Mark Jones, Samuel Martin, and Vincent Moulton. Below we outline the files in the different subfolders

## HIV
The HIV-1 sequence alignment is in the `hiv.fasta` file, originally from [1]. The Python script to construct the network in the Squirrel paper is in `hiv.py`, with the resulting network (in `eNewick` format) in the file `hiv.txt`.

## Primates
### All primates
The primate sequence alignment was too large for GitHub but can be retrieved as the file `1730_ALIGNMENT_CONCAT.paup.nex` from the [supplementary files](https://doi.org/10.5061/dryad.rfj6q577d) of the original source [2]. The Python script to construct the network in the Squirrel paper is in `primates.py`, with the resulting network (in `eNewick` format) in the file `primates.txt`.

### Cercopithecinae subfamily
The file `cerco.py` contains the script to construct the networks for the Cercopithecinae subfamily (using the same sequence alignment as the full primate data set). The two resulting networks in the Squirrel paper are (in `eNewick` format) in the file `cerco.txt`. The file `cerco_QNRSVM.csv` is from the paper [3] and contains the quarnets constructed with QNR-SVM (in their format). The file `cerco_QNRSVM.py` contains the script to generate the networks from their quarnets. 

## Xiphophorus
The tf-quarnets for the Xiphophorus alignment are in the zip-file `xiphophorus_quarnet_data.zip`. The tf-quarnets were generated with the MML-algorithm [4] and also come from that paper (using their format). The Python script to construct the network in the Squirrel paper is in `xiphophorus.py`, with the resulting network (in `eNewick` format) in the file `xiphophorus.txt`.

## Simulations
The networks for the simulation experiments in the Squirrel paper were generated with the script `generated_networks.py`. The generated networks are in the zip-file `generated_networks.zip` (using the `physquirrel` format for semi-directed networks). The filenames for the networks have the format `SDNXXXXX_LYYYYY_RZZZZZ.txt`, where `XXXXX` is the id-number, `YYYYY` the number of leaves, and `ZZZZZ` the number of reticulations. The numerical results of the experiments are in the file `experiment.csv` and are computed with the script `experiment.py`.

## References

1. M. Salemi and A.-M. Vandamme. The phylogenetic handbook: a practical approach to DNA and protein phylogeny. Cambridge University Press, 2003.
2. D. Vanderpool et al. “Primate phylogenomics uncovers multiple rapid radiations and ancient interspecific introgression”. In: PLOS Biology 18.12 (Dec. 2020), pp. 1–27. DOI: [10.1371/journal.pbio.3000954](https://doi.org/10.1371/journal.pbio.3000954).
3. T. Barton, E. Gross, C. Long, and J. Rusinko. Statistical learning with phylogenetic network invariants. 2022. arXiv:[2211.11919](https://arxiv.org/abs/2211.11919).
4. S. Martin, V. Moulton, and R. M. Leggett. Algebraic invariants for inferring 4-leaf semi-directed phylogenetic networks. 2023. bioRxiv:[10.1101/2023.09.11.557152](https://www.biorxiv.org/content/10.1101/2023.09.11.557152v3).
