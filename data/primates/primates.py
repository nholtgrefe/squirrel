from physquirrel import MSA

file = "/path/to/folder/1730_ALIGNMENT_CONCAT.paup.nex"
msa = MSA(data_folder + file)

# Make it a dense quarnet set
Q = msa.delta_heuristic(lam=0.3)

# Construct the network
outgroup = 'Mus_musculus'
networks, scores = Q.squirrel(include_score=True, all_networks=True, outgroup=outgroup)
networks[0].visualize(layout="dot",title="Constructed Primate network")
