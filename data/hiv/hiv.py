from physquirrel import MSA

# Load the alignment
file = '/path/to/folder/hiv.fasta'
msa = MSA(file)

# Create tf-quarnets
Q = msa.delta_heuristic(lam=0.3)

# Construct networks
networks, scores = Q.squirrel(include_score=True, triangles=False, all_networks=True, outgroup='C')
networks[0].visualize(layout='dot', title="Constructed HIV network")
