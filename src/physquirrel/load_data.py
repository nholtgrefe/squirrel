# Module containing functions to load quarnet data from different filetypes

import os, re
from .quarnetset import QuarnetSet
from .quarnet import SplitQuarnet, FourCycle, DoubleTriangle, SingleTriangle, QuartetTree
from .utils import CircularOrdering
from .splits import Split

def load_quarnets_from_invariants(folder, weight=True):
    """Returns a DenseQuarnetSet containing quarnets from a folder with invariant-files. 
    Optionally also extract weights from the data.
    Raises an error if a file is not in the correct format."""
    
    Q = QuarnetSet()
    for filename in os.listdir(folder):
        with open(folder + filename, 'r') as file:
            lines = file.readlines()

            tree_pattern = r'\(\(([^,]+),([^,]+)\),\(([^,]+),([^,]+)\)\)'
            tree_match = re.search(tree_pattern, lines[13])
            if tree_match:
                tree_lst = [tree_match.group(1), tree_match.group(2), tree_match.group(3), tree_match.group(4)]
                split = Split({tree_lst[0], tree_lst[1]}, {tree_lst[2],tree_lst[3]})
                quarnet = SplitQuarnet(split)
            else:
                cycle_pattern = r'\(([^,]+),([^,]+),([^,]+),([^,]+)\)'
                cycle_match = re.search(cycle_pattern, lines[1])
                if cycle_match:
                    cycle_lst = [cycle_match.group(1), cycle_match.group(2), cycle_match.group(3), cycle_match.group(4)]
                    circular_order = CircularOrdering(cycle_lst)
                    reticulation = cycle_lst[0]
                    quarnet = FourCycle(circular_order, reticulation)
                    
                    if weight:
                        weight_pattern = re.compile(r'\d+\.\d+e[+-]\d+')
                        values = []
                        for line in lines:
                            weight_match = weight_pattern.findall(line)
                            if len(weight_match) == 1:
                                values.append(float(weight_match[0]))
                        val = values[1] / values[0] - 1.0
                        quarnet.set_weight(min(val, 1.0))
                else:
                    raise ValueError("Could not extract quarnet from file.")
        Q.add_quarnet(quarnet)
        
    return Q.to_densequarnetset()




def load_quarnets_from_SVM(file, weight=True):
    """Returns a DenseQuarnetSet containing quarnets from a .xlsx file with SVM
    quarnet data. Optionally extracts weight from the data."""
    
    import pandas as pd
        
    quarnet_map = {
        1: QuartetTree(Split({1,2},{3,4})),
        2: QuartetTree(Split({1,3},{2,4})),
        3: QuartetTree(Split({2,3},{1,4})),
        4: SingleTriangle(Split({1,2},{3,4}), {3,4}),
        5: SingleTriangle(Split({1,3},{2,4}), {2,4}),
        6: SingleTriangle(Split({1,4},{2,3}), {2,3}),
        7: SingleTriangle(Split({2,3},{1,4}), {1,4}),
        8: SingleTriangle(Split({2,4},{1,3}), {1,3}),
        9: SingleTriangle(Split({3,4},{1,2}), {1,2}),
        10: FourCycle(CircularOrdering([2,3,4,1]), 2),
        11: FourCycle(CircularOrdering([4,3,2,1]), 4),
        12: FourCycle(CircularOrdering([2,4,3,1]), 2),
        13: FourCycle(CircularOrdering([3,4,2,1]), 3),
        14: FourCycle(CircularOrdering([2,3,1,4]), 2),
        15: FourCycle(CircularOrdering([1,3,2,4]), 1),
        16: FourCycle(CircularOrdering([1,4,3,2]), 1),
        17: FourCycle(CircularOrdering([3,4,1,2]), 3),
        18: FourCycle(CircularOrdering([1,3,4,2]), 1),
        19: FourCycle(CircularOrdering([4,3,1,2]), 4),
        20: FourCycle(CircularOrdering([3,2,4,1]), 3),
        21: FourCycle(CircularOrdering([4,2,3,1]), 4),
        22: DoubleTriangle(Split({1,2},{3,4}), {1}, {3}),
        23: DoubleTriangle(Split({1,3},{2,4}), {1}, {2}),
        24: DoubleTriangle(Split({2,3},{1,4}), {2}, {1})
        }
    
    df = pd.read_csv(file, delimiter=',')
    labels_lists = df.iloc[:, 0:4].values.tolist()
    network_types = [elt[0] for elt in df.iloc[:, 4:5].values.tolist()]
    if weight:
        network_weights = [float(elt[0]/100) for elt in df.iloc[:, 5:6].values.tolist()]

    filtered_df = df[df.iloc[:, 8].notna() & df.iloc[:, 9].notna()]
    mapping = dict(zip(filtered_df.iloc[:, 8], filtered_df.iloc[:, 9]))
    mapping = {int(key): value for key, value in mapping.items()}
    
    Q = QuarnetSet()
    for i in range(len(labels_lists)):
        q = quarnet_map[network_types[i]]
        labelling = labels_lists[i]
        m = {i:labelling[i-1] for i in [1,2,3,4]}
        q = q.relabel(m)
        q = q.relabel(mapping)
        if weight:
            q.set_weight(network_weights[i])
        Q.add_quarnet(q)
    
    return Q.to_densequarnetset()

    