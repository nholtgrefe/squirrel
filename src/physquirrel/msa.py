from itertools import combinations
import numpy as np
from functools import lru_cache

from .splits import Split
from .utils import CircularOrdering
from .quarnet import SplitQuarnet, FourCycle, CycleQuarnet
from .quarnetset import DenseQuarnetSet

class MSA:
    """
    Class for Multiple Sequence Alignments. If a .fasta, .fa, .nexus or .nex file is used as argument, that is used
    for initialization. Otherwise takes a dictionairy of sequences where keys are sequence IDs and values are sequences.
        self.taxa = a set of all taxa
        self.sequences = a dictionary of sequences where keys are sequence IDs and values are sequences.
        self.sequences_arrays = a dictionary where the sequences are np.arrays instead of strings
    """
    
    def __init__(self, arg):
        if type(arg) is str:
            if arg.endswith(('.fa', '.fasta')):
                self.sequences = self._read_fasta(arg)
            elif arg.endswith(('.nex', '.nexus')):
                self.sequences = self._read_nexus(arg)
            else:
                raise ValueError("Wrong filetype.")
        elif type(arg) is dict:
            self.sequences = arg
            
        self.taxa = set(self.sequences.keys())
        self.taxa_order = list(self.taxa)

    def __repr__(self):
        msa_repr = ""
        for seq_id, seq in self.sequences.items():
            msa_repr += f"> {seq_id:<30} {seq[0:40]} ...\n"
        return f"MSA on {len(self.taxa)} taxa [" + "\n" + msa_repr + "]"
    
    def __len__(self):
        """Returns the number of taxa in the MSA."""
        return len(self.taxa)
    
    @lru_cache(maxsize=1)
    def seq_length(self):
        """Returns the sequence length of the first sequence in the dictionary."""
        return len(list(self.sequences.values())[0])
    
    def delta_heuristic(self, lam=0.3, weight=True):
        """Create a DenseQuarnetSet from the MSA using the delta-heuristic.
        If weight is True, the quarnets will be weighted."""
        
        distance = self._compute_distance_matrix()
        delta_sum = {taxum:0 for taxum in self.taxa}
        
        quarnets = []
        index_combinations = combinations(range(len(self.taxa)), 4)
        
        for comb in index_combinations:
            i, j, k, l = list(comb)
            a = self.taxa_order[i]; b = self.taxa_order[j]; c = self.taxa_order[k]; d = self.taxa_order[l]

            split_distances = {
                Split({a,b},{c,d}): distance[i][j]+distance[k][l], 
                Split({a,c},{b,d}): distance[i][k]+distance[j][l], 
                Split({a,d},{b,c}): distance[i][l]+distance[j][k]
                }
            # sort the split distances in increasing order
            sorted_splits = sorted(split_distances.items(), key=lambda item: item[1])
            
            # compute the delta value
            if sorted_splits[0][1] == sorted_splits[1][1] and sorted_splits[1][1] == sorted_splits[2][1]:
                delta = 0
            else:
                numerator = sorted_splits[2][1] - sorted_splits[1][1]
                denominator = sorted_splits[2][1] - sorted_splits[0][1]
                delta = numerator / denominator
            
            delta_sum[a] += delta; delta_sum[b] += delta; delta_sum[c] += delta; delta_sum[d] += delta
            
            if delta < lam:
                best_split = sorted_splits[0][0]
                q = SplitQuarnet(best_split)
            else:
                wrong_split = sorted_splits[-1][0]
                if wrong_split == Split({a,b}, {c,d}):
                    circ = CircularOrdering([a,c,b,d])
                elif wrong_split == Split({a,c}, {b,d}):
                    circ = CircularOrdering([a,b,c,d])
                elif wrong_split == Split({a,d}, {b,c}):
                    circ = CircularOrdering([a,b,d,c])
                
                q = CycleQuarnet(circ)
            
            if weight == True:
                if delta < lam:
                    w = abs(lam - delta) / lam
                else:
                    w = abs(lam-delta) / (1 - lam)
                q.set_weight(w)
                
            quarnets.append(q)
        # Assign reticulations
        # order taxa according to (mean) delta values
        reticulation_order = sorted(delta_sum, key=lambda x: delta_sum[x], reverse=True)

        ret_quarnets = []
        for q in quarnets:
            if isinstance(q, CycleQuarnet):
                circ = q.circular_order
                ret = next((element for element in reticulation_order if element in q.leaves), None)
                
                if weight == True:
                    w = q.weight
                    q = FourCycle(circ, ret, weight=w)
                else:
                    q = FourCycle(circ, ret)

            ret_quarnets.append(q)          

        return DenseQuarnetSet(ret_quarnets)
    
    @lru_cache(maxsize=1)
    def _compute_distance_matrix(self):
        """Computes the matrix of distance-values for each pair of taxa, using the normalized Hamming-distance.
        Columns that contain a gap for eather of the two taxa are not considered. The ordering is
        as in self.taxa_order.
        See: delta Plots: A Tool for Analyzing Phylogenetic Distance Data
                by B.R. Holland, K.T. Huber, A. Dress, V. Moulton."""
        
        distance = np.zeros((len(self.taxa), len(self.taxa)))
        
        for i, a in enumerate(self.taxa_order):
            for j, b in enumerate(self.taxa_order):
                if i > j:
                    hamming = 0
                    no_gap_length = 0
                    
                    string_a = self.sequences[a]
                    string_b = self.sequences[b]
                    
                    for base1, base2 in zip(string_a, string_b):
                        allowed_charachters = ['A','C','G','T','a','c','t','g']
                        if base1 in allowed_charachters and base2 in allowed_charachters:
                            no_gap_length += 1
                            if base1.capitalize() != base2.capitalize():
                                hamming += 1

                    distance[i,j] = hamming / no_gap_length
                    distance[j,i] = hamming / no_gap_length
                
        
        return distance
    
    
    def _read_fasta(self, file_path):
        """Reads a fasta/fa file and returns dictionary of sequences."""

        sequences = {}
        with open(file_path, 'r') as file:
            sequence_id = None
            sequence_data = []

            for line in file:
                line = line.strip()

                if line.startswith('>'):
                    if sequence_id:
                        # Store the previous sequence
                        sequences[sequence_id] = ''.join(sequence_data)
                    
                    # Start a new sequence
                    sequence_id = line[1:]  # Remove the '>' character
                    sequence_data = []
                else:
                    sequence_data.append(line)

            # Store the last sequence
            if sequence_id:
                sequences[sequence_id] = ''.join(sequence_data)
        
        return sequences
    
    def _read_nexus(self, file_path):
        """Reads a nexus/nex file and returns dictionairy of sequences."""

        sequences = dict()
        with open(file_path, 'r') as file:
            in_matrix = False
            for line in file:
                line = line.strip()
                if line.lower().startswith('matrix'):
                    in_matrix = True
                    continue
                if in_matrix:
                    if line.lower().startswith(';'):
                        break
                    if line:
                        parts = line.split()
                        seq_name = parts[0]
                        seq_data = ''.join(parts[1:])
                        if seq_name in sequences:
                            sequences[seq_name] += seq_data
                        else:
                            sequences[seq_name] = seq_data
        return sequences      
