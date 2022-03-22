# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 11:53:32 2021

@author: 00ery
"""

#-----THIS CODE WORKS FOR DUMP FILES WITH 'mol x y z'-----

#importing modules
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

#increasing pixel density (300 should give very good quality)
sns.set(rc={"figure.dpi":100, 'savefig.dpi':100})


#constants
NO_OF_MOLECULES = 1024
NO_OF_BONDS = 29

BOND_DISTANCE_NC1 = 1.12
TOLERANCE_NC1 = 0.7  #0.7

BOND_DISTANCE_7S = 0.224
TOLERANCE_7S = 0.1

INDIRECT_BOND_DISTANCE_7S = 1
INDIRECT_TOLERANCE_7S = 0.4

BOX_SIZE = 32

#---FUNCTIONS---
#selecting all rows from dump file corresponding to given frame of simulation
def sim_frame(domain, file_name, n):
    if domain == "nc1":
        k = 1
    elif domain == "7s":
        k = 4
    elif domain == "all":
        k = 43
    elif domain == "bonds":
        k = NO_OF_BONDS
    #reading every line in file for 7s and nc1 cases
    file = open(file_name, 'r')
    file_lines = file.readlines()
    #selecting rows for given frame
    current_frame = file_lines[(9+(NO_OF_MOLECULES*k))*(n-1):(9+(NO_OF_MOLECULES*k))*(n)]
    current_frame = current_frame[9:len(current_frame)]
    file.close()
    
    return current_frame

#converting lines from file to position vectors (for nc1 and s7 only)
def pos(line):
    vector = np.array(line.split()[1:4])
    return np.array([float(x) for x in vector])

#converting lines from file to position vectors (for files with 'id type mol xyz')
def slice_line(line, start, end):
    vector = np.array(line.split()[start:end])
    return np.array([float(x) for x in vector])

#converting molecule index to integer (for nc1 and 7s only)
def mol(line):
    index = int(line.split()[0])
    return index


#calculating mean stress and orientation for given molecule (for files with 'id type mol c_stress[3] xyz)
def stress_and_orientation(frame):
    #first putting indices of main chain particles into list
    main_chain_indices_0 = list(range(1,30))
    main_chain_indices_all = []
    for i in range(1,NO_OF_MOLECULES+1):
        for j in main_chain_indices_0:
            main_chain_indices_all.append(j+((i-1)*43))
    
    #filtering out unneeded particles
    frame_main_chain = [i.split() for i in frame if int(slice_line(i, 0, 1)[0]) in main_chain_indices_all]
    
    #sorting by particle index
    frame_main_chain = sorted(frame_main_chain, key=lambda x: int(x[0]))
    
    #calculating mean orientation vector and stress for each molecule
    orientations = []
    stresses = []
    for i in range(NO_OF_MOLECULES):
        coords_single = []
        stresses_single = []
        for ii in frame_main_chain[i*29:(i+1)*29]:
            #retrieving stress and positions
            coords_single.append(np.array([float(x) for x in ii[4:7]]))
            stresses_single.append(float(ii[3]))
            #finding orientation vectors
            orientations_single = []
            for j in range(len(coords_single)-1):
                orientations_single.append(coords_single[j]-coords_single[j+1])
            #removing too long vectors which results from PBC
            orientations_single = [x for x in orientations_single if np.linalg.norm(x) < 20]
        orientations.append(np.mean(orientations_single, axis=0))
        stresses.append(np.mean(stresses_single))
    
    #normalising vectors
    orientations = [x/np.linalg.norm(x) for x in orientations]
    
    #quantifying alignment
    z_vector = np.array([0,0,1])
    alignments = [np.dot(z_vector, x)*np.dot(z_vector, x) for x in orientations]
          
    
    return alignments, stresses  

'''
#data = stress_and_orientation(sim_frame("all", "collagen_4_dump_SEED1000_CF1.0.xyz", 134))
sns.scatterplot([1,2,3,4,5], [1,2,3,4,5], label = "one", s=10)
sns.scatterplot([2,4,6,8,10], [3,5,7,9,11], label = "two", s=10)
plt.xlabel("$cos^2\u03B8$")
plt.ylabel("$T_{zz}$")
plt.legend()
'''

#finding fraction of binding sites which have formed a bond
def bonded_fraction(domain, frame):
    if domain == "nc1":
        k = 1
        dist = BOND_DISTANCE_NC1
        tol = TOLERANCE_NC1
    elif domain == "7s":
        k = 4
        dist = BOND_DISTANCE_7S
        tol = TOLERANCE_7S
    
    #saving only the coordinates of binding sites
    BS_pos = []
    for i in frame:
        BS_pos.append(pos(i))
    #counting the number of sites which are bonded to another site
    bonded_sites = 0
    for i in BS_pos:
        for ii in BS_pos:
            if abs((np.dot(i-ii,i-ii)**0.5)-dist) < tol:
                bonded_sites += 1
    return bonded_sites/(NO_OF_MOLECULES*k)


        
        
#finding the adjacency matrix
def adjacency_matrix(frame_nc1, frame_7s, matrix):
    #---NC1 domains---
    #saving binding site coordinates and molecule IDs
    BS_NC1_pos = []
    BS_NC1_mol = []
    for i in frame_nc1:
        BS_NC1_pos.append(pos(i))
        BS_NC1_mol.append(mol(i))
    #recording bonds in matrix    
    for i, j in zip(BS_NC1_pos, BS_NC1_mol):
        for ii, jj in zip(BS_NC1_pos, BS_NC1_mol):
            if abs((np.dot(i-ii,i-ii)**0.5)-BOND_DISTANCE_NC1) < TOLERANCE_NC1:
                matrix[j-1][jj-1] = 1
                matrix[jj-1][j-1] = 1
    
    #---7S domains---
    #saving binding site coordinates and molecule IDs
    BS_7S_pos = []
    BS_7S_mol = []
    for i in frame_7s:
        BS_7S_pos.append(pos(i))
        BS_7S_mol.append(mol(i))
    #recording bonds in matrix
    #directly bonded domains
    for i, j in zip(BS_7S_pos, BS_7S_mol):
        for ii, jj in zip(BS_7S_pos, BS_7S_mol):
            if abs((np.dot(i-ii,i-ii)**0.5)-BOND_DISTANCE_7S) < TOLERANCE_7S:
                matrix[j-1][jj-1] = 1
                matrix[jj-1][j-1] = 1
    #indirectly bonded domains
    for i, j in zip(BS_7S_pos, BS_7S_mol):
        for ii, jj in zip(BS_7S_pos, BS_7S_mol):
            if abs((np.dot(i-ii,i-ii)**0.5)-INDIRECT_BOND_DISTANCE_7S) < INDIRECT_TOLERANCE_7S and j!=jj:
                matrix[j-1][jj-1] = 1
                matrix[jj-1][j-1] = 1
    return matrix

#finding degree per node
def degree_per_node(G):
    edges = nx.number_of_edges(G)
    nodes = nx.number_of_nodes(G)
    return 2*edges/nodes

#finding mean betweenness centrality
def centrality(G):
    centrality_dict = nx.betweenness_centrality(G)
    centralities = []
    for key in centrality_dict:
        centralities.append(centrality_dict[key])
    return np.mean(centralities)

#mean cycle length
def cycle_length(G):
    total = 0
    for i in range(NO_OF_MOLECULES):
        try:
            total += len(nx.find_cycle(G,i))
        except:
            total += 0
    return total/NO_OF_MOLECULES

#writing contents of a list to a file
def list_to_file(file_name, list_name):
    #creating file
    f = open(file_name, "w")
    #function for writing lines to file
    def wl(text):
        f.write(str(text)+"\n")
    #writing each item in list to the file
    for i in list_name:
        wl(i)
    f.close()

#writing contents of file to list (for floating point numbers)  
def file_to_list(file_name):
    f = open(file_name, 'r')
    f_lines = f.readlines()
    for i in range(len(f_lines)):
        f_lines[i] = float(f_lines[i])
    return f_lines
  
#converting to float and taking mean
def mean_force(frame):
    float_frame = [float(i) for i in frame]
    return np.mean(float_frame)

#wiener index of largest structure of a graph
def wiener_index(G):
    #finding largest structure
    largest_cc = max(nx.connected_components(G), key=len)
    S = G.subgraph(largest_cc).copy()
    #returning wiener index
    return nx.wiener_index(S)

#outputting array with bonded/unbonded status of nc1 domain
def adjacency_array_nc1(frame):
    #initializing array
    nc1_array = np.zeros((NO_OF_MOLECULES))
    #saving binding site coordinates and molecule IDs
    BS_NC1_pos = []
    BS_NC1_mol = []
    for i in frame:
        BS_NC1_pos.append(pos(i))
        BS_NC1_mol.append(mol(i))
    #recording bonds in array    
    for i, j in zip(BS_NC1_pos, BS_NC1_mol):
        for ii, jj in zip(BS_NC1_pos, BS_NC1_mol):
            if abs((np.dot(i-ii,i-ii)**0.5)-BOND_DISTANCE_NC1) < TOLERANCE_NC1:
                nc1_array[j-1] = 1
                nc1_array[jj-1] = 1
    return nc1_array
    
#returning number of broken and formed bonds between a pair of adjacency arrays
def formed_broken_nc1(array_t, array_t_minus_1):
    delta_array = array_t - array_t_minus_1
    formed_bonds = ((delta_array == 1).sum())
    broken_bonds = ((delta_array == -1).sum())
    return [formed_bonds, broken_bonds]
 
    


'''
delta_bonds = []
broken_bonds = []
formed_bonds = []
total_bonds = []
for t in range(120,140,1):
    frame1 = sim_frame("nc1", "collagen_4_NC1_dump_SEED1000_CF1.0.xyz", t)
    ar1 = adjacency_array_nc1(frame1)
    frame2 = sim_frame("nc1", "collagen_4_NC1_dump_SEED1000_CF1.0.xyz", t+1)
    ar2 = adjacency_array_nc1(frame2)
    delta_bonds.append(formed_broken_nc1(ar2,ar1)[0]-formed_broken_nc1(ar2,ar1)[1])
    broken_bonds.append(formed_broken_nc1(ar2,ar1)[1])
    formed_bonds.append(formed_broken_nc1(ar2,ar1)[0])
    total_bonds.append(np.count_nonzero(ar1))
    print(t)
    
plt.plot(total_bonds)

total_bonds2 = [total_bonds[0]]
for i in delta_bonds:
    total_bonds2.append(total_bonds2[-1]+i)
    
plt.plot(total_bonds2)
'''


''' 
op = []
op_var = []  
for i in range(1,715,10):
    data = stress_and_orientation(sim_frame("all", "collagen_4_dump_SEED1000_CF1.0.xyz", i))[0]
    op.append(np.mean(data))
    op_var.append(np.std(data))
    print(i)
list_to_file("OP_STDEV.txt", op_var)
list_to_file("OP.txt", op)    
'''

'''
#PROPERTY AS FUNCTION OF TIME
no_of_frames = 100
no_of_seeds = 10

quantity1 = np.zeros((no_of_frames))
quantity1_STDEV = np.zeros((no_of_frames))

quantity2 = np.zeros((no_of_frames))
quantity2_STDEV = np.zeros((no_of_frames))

for k in range(100):
    total_quantity1 = np.zeros((no_of_seeds))
    total_quantity2 = np.zeros((no_of_seeds))
    for i in range(no_of_seeds):
        
        frame_t = sim_frame("nc1", f"collagen_4_NC1_dump_SEED{1000+i}_CF1.0.xyz", k+100)
        frame_t_minus_1 = sim_frame("nc1", f"collagen_4_NC1_dump_SEED{1000+i}_CF1.0.xyz", k-1+100)
        
        total_quantity1[i] = formed_broken_nc1(adjacency_array_nc1(frame_t),adjacency_array_nc1(frame_t_minus_1))[0]
        total_quantity2[i] = formed_broken_nc1(adjacency_array_nc1(frame_t),adjacency_array_nc1(frame_t_minus_1))[1]
        
    quantity1[k] = np.mean(total_quantity1)
    quantity1_STDEV[k] = np.std(total_quantity1)
    
    quantity2[k] = np.mean(total_quantity2)
    quantity2_STDEV[k] = np.std(total_quantity2)
    
    print(f"Current frame: {k}")
    

 
list_to_file("formed100200_STDEV.txt", quantity1_STDEV)
list_to_file("formed100200.txt", quantity1)
list_to_file("broken100200_STDEV.txt", quantity2_STDEV)
list_to_file("broken100200.txt", quantity2)
'''

    




'''
#DRAWING GRAPH REPRESENTATION OF NETWORKS
frame_nc1 = sim_frame(f"nc1", f"collagen_4_NC1_dump_SEED1009_CF1.0.xyz", 140)
frame_7s = sim_frame(f"7s", f"collagen_4_7S_dump_SEED1009_CF1.0.xyz", 140)
A = np.zeros((NO_OF_MOLECULES,NO_OF_MOLECULES))
A = adjacency_matrix(frame_nc1, frame_7s, A)
G = nx.from_numpy_matrix(A)
nx.draw_kamada_kawai(G, node_size=5)
'''




'''
#FORCE ON BONDS AS FUNCTION OF TIME
stress = np.zeros((463))
stress_STDEV = np.zeros((463))

for k in range(462):
    total_stress = np.zeros((10))
    
    for i in range(10):
        total_stress[i] = mean_force(sim_frame(f"bonds", f"force_SEED{1000+i}_CF1.0.xyz", k))
    
    stress[k] = np.mean(total_stress)
    stress_STDEV[k] = np.std(total_stress)
    
    print(f"Current frame: {k}")


list_to_file("force_SRshort.txt", stress)
list_to_file("force_SRshort_STDEV.txt", stress_STDEV)
'''






'''
#GRAPH PROPERTY AS FUNCTION OF CONCENTRATION
#selecting frame of interest (most likely the final one)
k = 120

qnt = np.zeros((10))
qnt_STDEV = np.zeros((10))

#iterating over concentrations
for cf in range(100,200,10):
    total_qnt = np.zeros((5))
    #iterating over seed
    for i in range(5): 
        A = np.zeros((NO_OF_MOLECULES,NO_OF_MOLECULES))
        A = adjacency_matrix(sim_frame("nc1", f"collagen_4_NC1_dump_SEED{1000+i}_CF{cf/100}.xyz", k), sim_frame("7s", f"collagen_4_7S_dump_SEED{1000+i}_CF{cf/100}.xyz", k), A)
        G = nx.from_numpy_matrix(A)
        total_qnt[i] = wiener_index(G)
    qnt[range(100,200,10).index(cf)] =  np.mean(total_qnt)
    qnt_STDEV[range(100,200,10).index(cf)] =  np.std(total_qnt)
    print(cf/100)
    
 
#writing quantity as function of time (for given concentration) to a file
list_to_file(f"wiener_of_CONC.txt", qnt)
list_to_file(f"wiener_of_CONC_STD.txt", qnt_STDEV)
'''


'''
#GRAPH PROPERTY AS FUNCTION OF ENERGY
#selecting frame of interest (most likely the final one)
k = 120

qnt = np.zeros((12))
qnt_STDEV = np.zeros((12))
energies = range(2,26,2)
#iterating over concentrations
for energy in energies:
    total_qnt = np.zeros((5))
    #iterating over seed
    for i in range(5): 
        A = np.zeros((NO_OF_MOLECULES,NO_OF_MOLECULES))
        A = adjacency_matrix(sim_frame("nc1", f"collagen_4_NC1_dump_SEED{1000+i}_CF1.0_E{energy}.xyz", k), sim_frame("7s", f"collagen_4_7S_dump_SEED{1000+i}_CF1.0_E{energy}.xyz", k), A)
        G = nx.from_numpy_matrix(A)
        total_qnt[i] = wiener_index(G)
    qnt[energies.index(energy)] =  np.mean(total_qnt)
    qnt_STDEV[energies.index(energy)] =  np.std(total_qnt)
    print(energy)
    
 
#writing quantity as function of time (for given concentration) to a file
list_to_file(f"wiener_of_E.txt", qnt)
list_to_file(f"wiener_of_E_STDEV.txt", qnt_STDEV)

'''




data = stress_and_orientation(sim_frame("all", "collagen_4_dump_SEED1000_CF1.0.xyz", 1000))[0]
print(np.mean(data))