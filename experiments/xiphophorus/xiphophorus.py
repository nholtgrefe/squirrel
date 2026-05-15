from physquirrel import DenseQuarnetSet, load_quarnets_from_invariants
import random
random.seed(42)

# Folder with quarnet data from MML algorithm
data_folder = "/path/to/folder/xiphophorus_quarnet_data/

# Load the quarnet set
Q = load_quarnets_from_invariants(data_folder, weight=True)

# We remove the taxon 'Xbirref' since it is not a real species an an artifact of the alignment procedure
# We remove one of the outgroups since Squirrel can only use one outgroup to root
keep = {#'Priapella',
 'Psjonesii',
 'Xalvarezi',
 'Xandersi',
 'Xbirchmann',
 #'Xbirref',
 'Xclemencia',
 'Xcontinens',
 'Xcortezi',
 'Xcouchianu',
 'Xevelynae',
 'Xgordoni',
 'Xhellerii',
 'Xmaculatus',
 'Xmalinche_',
 'Xmayae',
 'Xmeyeri',
 'Xmilleri',
 'Xmontezuma',
 'Xmonticolu',
 'Xmultiline',
 'Xnezahuaco',
 'Xnigrensis',
 'Xpygmaeus',
 'Xsignum',
 'Xvariatus',
 'Xxiphidium'}

newQ = set()
for q in Q:
    if len(q.leaves & keep) == 4:
        newQ.add(q)

# Make it a dense quarnet set
DQ = DenseQuarnetSet(newQ)

# Construct the networks
outgroup = 'Psjonesii'
networks, scores = DQ.squirrel(include_score=True, all_networks=True, outgroup=outgroup)
networks[0].visualize(layout="dot",title="Constructed Xiphophorus network")
