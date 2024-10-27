from physquirrel import MSA

# Load sequence alignment
file = "/path/to/folder/1730_ALIGNMENT_CONCAT.paup.nex"
msa = MSA(file)

# Restrict sequence alignment to Cercopithecinae subfamily (with outgroup Colobus angolensis pallatus)
seq = msa.sequences
keep = {
  'Cercocebus_atys',
  'Chlorocebus_sabaeus',
  'Colobus_angolensis_palliatus',
  'Macaca_fascicularis',
  'Macaca_mulatta',
  'Macaca_nemestrina',
  'Mandrillus_leucophaeus',
  'Papio_anubis',
  'Theropithecus_gelada',
 }

new_seq = {k:v for k,v in seq.items() if k in keep}
msa = MSA(new_seq)

# Create tf-quarnets
Q = msa.delta_heuristic(lam=0.3)

# Construct networks
outgroup = 'Colobus_angolensis_palliatus'
networks, scores = Q.squirrel(include_score=True, verbose=True, triangles=False, all_networks=True, outgroup=outgroup)

networks[0].visualize(layout="dot",title="Constructed Cercopithecinae network 1")
networks[1].visualize(layout="dot",title="Constructed Cercopithecinae network 2")
