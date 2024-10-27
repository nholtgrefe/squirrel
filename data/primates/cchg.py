from physquirrel import MSA

data_folder = "/home/nholtgreve/Documents/Projects/Phylogenetics/Level1 heuristic from quarnets/Code/SQuIrReL/data/primate/"
file ='primate.nex'
file = ""

msa = MSA(data_folder + file)

seq = msa.sequences

# Remove one of the outgroups (our algorithm uses a single outgroup to root)
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

# Make it a dense quarnet set
DQ = msa.delta_heuristic(lam=0.3)

# Construct the network
outgroup = 'Colobus_angolensis_palliatus'
networks, scores = DQ.squirrel(include_score=True, verbose=True, triangles=False, all_networks=True, outgroup=outgroup)
networks[0].visualize(layout="dot",title="Constructed Primate network 1")

networks[1].visualize(layout="dot",title="Constructed Primate network 2")
