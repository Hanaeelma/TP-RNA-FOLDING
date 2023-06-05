def load():
    print('src file loaded !')

import pandas as pd

def read_pdb_file(filename):
    pdb_dict = {}
    i = 0
    with open(filename, 'r') as pdb_file:
        for line in pdb_file:
            i += 1
            if line.startswith('ATOM'):
                record_type = line[0:6].strip()
                atom_number = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip()
                chain_id = line[21:22].strip()
                residue_number = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                
                if i not in pdb_dict:
                    if 'C3' in  atom_name:
                        pdb_dict[i] = {'record_type': record_type,
                            'atom_number': atom_number,
                            'atom_name': atom_name,
                            'residue_name': residue_name,
                            'residue_number': residue_number,
                            'x': x,
                            'y': y,
                            'z': z,
                        }
                
    return pdb_dict

import math

def calculate_distance(atom1, atom2): 
    x1, y1, z1 = atom1['x'], atom1['y'], atom1['z']
    x2, y2, z2 = atom2['x'], atom2['y'], atom2['z']
    distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
    return distance

def calculate_ca_distances(pdb_dict):
    distances = {}
    for key1, atom1 in pdb_dict.items():
        for key2, atom2 in pdb_dict.items():
            distance = calculate_distance(atom1, atom2)
            couple = atom1['residue_name'] + '-' + atom2['residue_name'] 
            if couple not in distances:
                distances[couple] = []
            if distance<20: #filtring big distances
                distances[couple].append(distance)
    return distances



## Step 2
import pandas as pd
import numpy as np


def count_occurrences(distance, residue_pairs, intervals):
    N = pd.DataFrame(index=intervals, columns=residue_pairs).fillna(0)

    for interval in intervals:
        for residue_pair in residue_pairs:
            if residue_pair in distance:
                for v in distance[residue_pair]:
                    if interval <= v <= interval + 1:
                        N.loc[interval, residue_pair] += 1

    return N


def calculate_Nxx(N):
    return N.sum(axis=0)


def calculate_Nij(N):
    return N.sum(axis=1)


def calculate_observed_frequencies(N, Nij):
    return N.divide(Nij, axis=0)


def calculate_reference_frequencies(N, Nxx):
    return N.divide(Nxx, axis=1)


def calculate_scores(obs_freq, ref_freq):
    score = -1 * np.log10(obs_freq / ref_freq)
    score = score.clip(upper=10).fillna(10)
    return score



import pandas as pd
import numpy as np

def calculate_frequencies(distance):

    residue_pairs = ['A-A', 'A-U', 'A-C', 'A-G', 'U-U', 'U-G', 'U-C', 'C-C', 'C-G', 'G-G']
    intervals = range(0, 21)
    
    N = pd.DataFrame(index=intervals, columns=residue_pairs).fillna(0)

    for interval in intervals:
        for residue_pair in residue_pairs:
            if residue_pair in distance:
                for v in distance[residue_pair]:
                    if interval <= v <= interval + 1:
                        N.loc[interval, residue_pair] += 1

    Nxx = N.sum(axis=0)
    Nij = N.sum(axis=1)
    
    obs_freq = N.divide(Nij, axis=0)
    ref_freq = N.divide(Nxx, axis=1)

    score = -1 * np.log10(obs_freq / ref_freq)
    score = score.clip(upper=10).fillna(10)


    for s in score.columns:
        outputname = './data/scores/'+s+'.csv'
        score[s].to_csv(outputname , index=False, header=True)
    
    return obs_freq, ref_freq, score
