from physquirrel import SemiDirectedNetwork
import os, time

# Epsilon value to perturb quarnets
epsilons = [0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]

# Folder containing generated networks
data_folder = "/home/nholtgreve/Documents/Projects/Phylogenetics/Level1 heuristic from quarnets/Code/Lev1TSP/data/simulations/240927/"
files = os.listdir(data_folder)
files.sort()

result_file = "/home/nholtgreve/Documents/Projects/Phylogenetics/Level1 heuristic from quarnets/Code/Lev1TSP/data/simulations/result_experiment.txt"

with open(result_file, 'a') as f:
    f.write("file \t epsilon \t QC \t QS \t QCprime \t ret \t time \n")

for file in files:
    N = phylo.SemiDirectedNetwork()
    N.load_from_file(data_folder + file)
    Q = N.quarnets(triangles=False)
    
    for epsilon in epsilons:
        Qpert = Q._shake(epsilon)
        
        time1 = time.perf_counter()
        N_reconstructed = Qpert.reconstruct_network(verbose=False, measure="QC", triangles=False)
        time2 = time.perf_counter()
        
        Q_recon = N_reconstructed.quarnets(triangles=False)
        
        score1 = Q.consistency(Q_recon, triangles=False) # QC
        score2 = Q.distance(Q_recon, triangles=False) #QS
        score3 = Q.consistency_prime(Q_recon, reference=Qpert, triangles=False) #QCprime
        
        ret = len(N_reconstructed.reticulation_nodes())

        timing = time2 - time1

        info = [epsilon, score1, score2, score3, ret, timing]
        with open(result_file, 'a') as f:
            f.write(file + "\t")
            f.write("\t".join(map(str, info)) + "\n")
        print(file, epsilon)
