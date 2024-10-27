from physquirrel import load_quarnets_from_SVM

# Load the quarnets
file = '/path/to/folder/cerco_QNRSVM.csv
DQ = load_quarnets_from_SVM(file)

# Adjust weights for well-supported quarnets as in the paper
w = 0.95
for q in DQ:
    if q.weight < w:
        q.set_weight(0.0)
    elif q.weight >= w:
        q.set_weight(1.0)

# Construct the networks
outgroup = 'Colobus_angolensis_palliatus'
networks, scores = DQ.squirrel(include_score=True, triangles=True, all_networks=True, outgroup=outgroup)
networks[0].visualize(font_size=7, title="Constructed Cercopithecinae network")
print(scores[0])
