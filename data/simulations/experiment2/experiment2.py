import physquirrel as psq
import os, re

# Folder containing orginal semi-directed networks
network_folder = "/path/to/folder/with/generated_sdnetworks/"
network_files = os.listdir(network_folder)
network_files.sort()

# Folder containing simulated alignments
sequence_folder = "/path/to/folder/with/simulated_alignments/"
sequence_files = os.listdir(sequence_folder)
sequence_files.sort()

# File to write results to
result_file = "/path/to/file/with/results/results_experiment2.txt"

with open(result_file, 'a') as f:
    f.write("file \t QC \t QS \t ret \n")

for network_file in network_files:
    match = re.match(r"SDN(\d+)_L(\d+)_R(\d+)", network_file)
    id_string, _, _ = match.groups()

    # Load original SD-network
    N = psq.SemiDirectedNetwork()
    N.load_from_file(network_folder + network_file)

    # Create quarnets of original network
    Q = N.quarnets(triangles=False)

    # Get filenames of corresponding alignments
    sequence_files_id = [f for f in sequence_files if f.startswith('SEQ'+id_string)]

    for sequence_file in sequence_files_id:
        # Create an MSA object from alignment
        msa = psq.MSA(sequence_folder + sequence_file)

        # Use delta-heuristic to create quarnets
        Q_delta = msa.delta_heuristic()

        # Use Squirrel to construct network
        N_reconstructed = Q_delta.squirrel()

        # Get quarnets of constructed network
        Q_recon = N_reconstructed.quarnets()

        #Compute the two scores
        score1 = Q.consistency(Q_recon, triangles=False) # QC
        score2 = 1 - Q.distance(Q_recon, triangles=False) #QS

        # Get number of reticulations of newly constructed network
        ret = len(N_reconstructed.reticulation_nodes())

        # Write results to file
        info = [score1, score2, ret]
        with open(result_file, 'a') as f:
            f.write(sequence_file + "\t")
            f.write("\t".join(map(str, info)) + "\n")
