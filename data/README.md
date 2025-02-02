# Supplementary data and scripts
This folder contains the biological datasets, generated networks, simulated sequences, python scripts and results of the simulation experiments of the paper:
> **Squirrel: Reconstructing semi-directed phylogenetic level-1 networks from four-leaved networks or sequence alignments**.
> *Niels Holtgrefe, Katharina T. Huber, Leo van Iersel, Mark Jones, Samuel Martin, and Vincent Moulton.* 2024. bioRxiv: [10.1101/2024.11.01.621567](https://doi.org/10.1101/2024.11.01.621567)

The results in the paper and in this repository were all obtained using version `1.0.6` of `physquirrel`.

## HIV
The HIV-1 sequence alignment is in the `hiv.fasta` file, originally from [1]. The Python script to construct the network in the Squirrel paper is in `hiv.py`, with the resulting network (in `eNewick` format) in the file `hiv.txt`.

## Primates
### All primates
The primate sequence alignment is too large for GitHub but can be retrieved as the file `1730_ALIGNMENT_CONCAT.paup.nex` from the [supplementary files](https://doi.org/10.5061/dryad.rfj6q577d) of the original source [2]. The Python script to construct the network in the Squirrel paper is in `primates.py`, with the resulting network (in `eNewick` format) in the file `primates.txt`.

### Cercopithecinae subfamily
The file `cerco.py` contains the script to construct the networks for the Cercopithecinae subfamily (using the same sequence alignment as the full primate data set). The two resulting networks in the Squirrel paper are (in `eNewick` format) in the file `cerco.txt`. The file `cerco_QNRSVM.csv` is from the paper [3] and contains the quarnets constructed with QNR-SVM (in their format). The file `cerco_QNRSVM.py` contains the script to generate the networks from their quarnets. 

## Xiphophorus
The tf-quarnets for the Xiphophorus alignment are in the zip-file `xiphophorus_quarnet_data.zip`. The tf-quarnets were generated with the MML-algorithm [4] and also come from that paper (using their format). The Python script to construct the network in the Squirrel paper is in `xiphophorus.py`, with the resulting network (in `eNewick` format) in the file `xiphophorus.txt`.

## Simulations
### Experiment 1
The semi-directed networks for the first simulation experiment in the Squirrel paper were generated with the script `generated_networks.py`. The generated networks are in the zip-file `generated_sdnetworks.zip` (using the `physquirrel` format for semi-directed networks). The filenames for the networks have the format `SDNXXXXX_LYYYYY_RZZZZZ.txt`, where `XXXXX` is the id-number, `YYYYY` the number of leaves, and `ZZZZZ` the number of reticulations. The numerical results of the experiments are in the file `results_experiment1.csv` and are computed with the script `experiment1.py`.

### Experiment 2
The directed networks with branch lengths used to simulate sequences in the second simulation experiment of the Squirrel paper were generated with part of the script `simulate_alignments.py`. These generated directed networks are obtained from the semi-directed networks of the first experiment and are in the zip-file `generated_dnetworks.zip` (using the eNewick format). The filenames for the networks have the format `WDNXXXXX_LYYYYY_RZZZZZ.txt`, where `XXXXX` is the id-number the corresponding semi-directed network, `YYYYY` the number of leaves, and `ZZZZZ` the number of reticulations. 

The sequence alignments that were simulated for the second experiment are available .... (in `.nexus` format). The alignments were simulated using `SeqGen` [5] in the script `simulate_alignments.py`. The sequence alignments are in the format `SEQXXXXX_LYYYYY_RZZZZZ_LENAAAAAAAAAAAA.nexus`, where `XXXXX` is the id-number the corresponding (semi-)directed network, `YYYYY` the number of leaves, `ZZZZZ` the number of reticulations, and `AAAAAAAAAAAA` the length (in bp) of the sequences. The numerical results of the experiment are in the file `results_experiment2.csv` and are computed with the script `experiment2.py`.


## References

1. M. Salemi and A.-M. Vandamme. The phylogenetic handbook: a practical approach to DNA and protein phylogeny. Cambridge University Press, 2003.
2. D. Vanderpool et al. “Primate phylogenomics uncovers multiple rapid radiations and ancient interspecific introgression”. In: PLOS Biology 18.12 (Dec. 2020), pp. 1–27. DOI: [10.1371/journal.pbio.3000954](https://doi.org/10.1371/journal.pbio.3000954).
3. T. Barton, E. Gross, C. Long, and J. Rusinko. Statistical learning with phylogenetic network invariants. 2022. arXiv:[2211.11919](https://arxiv.org/abs/2211.11919).
4. S. Martin, V. Moulton, and R. M. Leggett. Algebraic invariants for inferring 4-leaf semi-directed phylogenetic networks. 2023. bioRxiv:[10.1101/2023.09.11.557152](https://www.biorxiv.org/content/10.1101/2023.09.11.557152v3).
5. A. Rambaut, N. C. Grass. Seq-gen: an application for the monte carlo simulation of dna sequence evolution along phylogenetic trees. 1997. Bioinformatics. 13:235–238. DOI: [10.1093/bioinformatics/13.3.235](https://doi.org/10.1093/bioinformatics/13.3.235).
