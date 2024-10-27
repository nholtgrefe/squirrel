from physquirrel import random_semi_directed_network
import random, math

def generate():
    data_folder = "/path/to/target/folder/to/save/networks/"
    
    leaf_number_list = [10,15,20,25,30,35]
    nr_networks = 100
    
    id_nr = 0
    for leaf_nr in leaf_number_list:
        for net_nr in range(nr_networks):
            
            max_ret = math.floor(math.floor(leaf_nr/3))
            ret_nr = random.randint(0, max_ret)
            
            while True:
                try:
                    N = random_semi_directed_network(leaf_nr, ret_nr)
                    break
                except:
                    pass
                
            id_string = str(id_nr+1).zfill(5)
            leaf_string = str(leaf_nr).zfill(5)
            ret_string = str(ret_nr).zfill(5)
    
            file = f"SDN{id_string}_L{leaf_string}_R{ret_string}.txt"
            N.save_to_file(data_folder + file)
    
            id_nr += 1
