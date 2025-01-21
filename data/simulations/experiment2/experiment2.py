import physquirrel as psq
import os, time, re
    
network_folder = "/home/nholtgreve/Documents/Projects/Phylogenetics/Level1 heuristic from quarnets/Code/SQuIrReL/data/experiments/generated_networks/"
network_files = os.listdir(network_folder)
network_files.sort()

sequence_folder = "/home/nholtgreve/Documents/Projects/Phylogenetics/Level1 heuristic from quarnets/Code/SQuIrReL/data/experiments/simulated_alignments/"
sequence_files = os.listdir(sequence_folder)
sequence_files.sort()

result_file = "/home/nholtgreve/Documents/Projects/Phylogenetics/Level1 heuristic from quarnets/Code/SQuIrReL/data/experiments/result_sequence_experiment.txt"

with open(result_file, 'a') as f:
    f.write("file \t delta-QC \t delta-1-QS \t squirrel-QC \t squirrel-1-QS \t ret \t delta-time \t squirrel-time \n")

for network_file in network_files:
    match = re.match(r"SDN(\d+)_L(\d+)_R(\d+)", network_file)
    id_string, _, _ = match.groups()
    
    N = psq.SemiDirectedNetwork()
    N.load_from_file(network_folder + network_file)

    Q = N.quarnets(triangles=False)
    
    sequence_files_id = [f for f in sequence_files if f.startswith('SEQ'+id_string)]

    for sequence_file in sequence_files_id:
                    
        msa = psq.MSA(sequence_folder + sequence_file)
        
        time1 = time.perf_counter()
        Q_delta = msa.delta_heuristic()
        time2 = time.perf_counter()
                    
        N_reconstructed = Q_delta.squirrel()
        time3 = time.perf_counter()
        
        Q_recon = N_reconstructed.quarnets()
        
        score1 = Q.consistency(Q_delta, triangles=False) # QC
        score2 = 1 - Q.distance(Q_delta, triangles=False) #QS
        
        score3 = Q.consistency(Q_recon, triangles=False) # QC
        score4 = 1 - Q.distance(Q_recon, triangles=False) #QS
        
        time_delta = time2 - time1
        time_squirrel = time3 - time2
        
        ret = len(N_reconstructed.reticulation_nodes())

        info = [score1, score2, score3, score4, ret, time_delta, time_squirrel]
        with open(result_file, 'a') as f:
            f.write(sequence_file + "\t")
            f.write("\t".join(map(str, info)) + "\n")
