This folder contains the biological datasets, the simulated networks and results of the simulation experiments from the paper:
'Squirrel: Reconstructing semi-directed phylogenetic level-1 networks from four-leaved networks and sequence alignments' by Niels Holtgrefe, Katharina T. Huber, Leo van Iersel, Mark Jones, Samuel Martin, and Vincent Moulton. Below we outline the files in the different subfolders

## HIV
The HIV-1 sequence alignment is in the `hiv.fasta` file, originally from [1]. The Python script to construct the network in the Squirrel paper is in `hiv.py`, with the resulting network (in `eNewick` format) in the file `hiv.txt`.

## Primates
### Full primate set
The primate sequence alignment was too large for GitHub but can be retrieved as the file `1730_ALIGNMENT_CONCAT.paup.nex` from the [supplementary files](https://doi.org/10.5061/dryad.rfj6q577d) of the original source [2]. The Python script to construct the network in the Squirrel paper is in `primates.py`, with the resulting network (in `eNewick` format) in the file `primates.txt`.


## Xiphophorus
The tf-quarnets for the Xiphophorus alignment are in the zip-file `xiphophorus_quarnet_data.zip`. The tf-quarnets were generated with the MML-algorithm [4] and also come from that paper (using their format). The Python script to construct the network in the Squirrel paper is in `xiphophorus.py`, with the resulting network (in `eNewick` format) in the file `xiphophorus.txt`.

## Simulations
The networks generated for the simulation experiments in the paper are in the file `dsfdsf`. The filenames have the format `oiij;oji`. The numerical results of the experiments are in the file `ddfgd`

## References

1. M. Salemi and A.-M. Vandamme. The phylogenetic handbook: a practical approach to DNA and protein phylogeny. Cambridge University Press, 2003.
2. D. Vanderpool et al. “Primate phylogenomics uncovers multiple rapid radiations and ancient interspecific introgression”. In: PLOS Biology 18.12 (Dec. 2020), pp. 1–27. DOI: 10.1371/journal.pbio.3000954.
3. T. Barton, E. Gross, C. Long, and J. Rusinko. Statistical learning with phylogenetic network invariants. 2022. arXiv:2211.11919.
4. S. Martin, V. Moulton, and R. M. Leggett. Algebraic invariants for inferring 4-leaf semi-directed phylogenetic networks. 2023. bioRxiv: 10.1101/2023.09.11.557152.
