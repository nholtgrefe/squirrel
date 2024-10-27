from physquirrel import SemiDirectedNetwork
import os, time

# Epsilon value to perturb quarnets
epsilons = [0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]

# Folder containing generated networks
data_folder = "/folder/containing/generated_networks/"
files = os.listdir(data_folder)
files.sort()

# File to save numerical results
result_file = "/path/to/file/with/result/result_experiment.txt"

# Create a header line
with open(result_file, 'a') as f:
    f.write("file \t epsilon \t QC \t 1-QS \t ret \t time \n")

# Iterate through networks
for file in files:
    N = SemiDirectedNetwork()
    N.load_from_file(data_folder + file)

    # Generate tf-quarnets of the network
    Q = N.quarnets(triangles=False)
    
    for epsilon in epsilons:
        # Perturb a fraction of epsilon of the tf-quarnets
        Qpert = Q._shake(epsilon)

        # Construct new network from perturbed quarnets
        time1 = time.perf_counter()
        N_reconstructed = Qpert.reconstruct_network(verbose=False, measure="QC", triangles=False)
        time2 = time.perf_counter()

        # Retrieve results
        Q_recon = N_reconstructed.quarnets(triangles=False)
        
        score1 = Q.consistency(Q_recon, triangles=False) # QC
        score2 = Q.distance(Q_recon, triangles=False) #QS
        
        ret = len(N_reconstructed.reticulation_nodes())

        timing = time2 - time1

        # Write results to file
        info = [epsilon, score1, score2, ret, timing]
        with open(result_file, 'a') as f:
            f.write(file + "\t")
            f.write("\t".join(map(str, info)) + "\n")
        print(file, epsilon)

# This code will produce a tab-delimited .txt file, which we later converted into a .csv file for better viewing on GitHub
