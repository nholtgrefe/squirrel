from physquirrel import MSA

# Load the alignment
file = "/path/to/folder/1730_ALIGNMENT_CONCAT.paup.nex"
msa = MSA(data_folder + file)

# Create tf-quarnets
Q = msa.delta_heuristic(lam=0.3)

# Construct networks
outgroup = 'Mus_musculus'
networks, scores = Q.squirrel(include_score=True, all_networks=True, outgroup=outgroup)
networks[0].visualize(layout="dot",title="Constructed Primate network")
