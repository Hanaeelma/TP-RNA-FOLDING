# Project_RNA_Prediction
# TP RNA
For a given ribonucleotide chain, the RNA folding problem consists in finding the native fold
among the astronomically large number of possible conformations. The native fold being the
one with the lowest Gibbs free energy, the objective function should be an estimator of this
energy.
```{python}
# install packages
! pip install -r requirements.txt
```
### Part 1
Train the objective function, using interatomic distance distributions that are computed from a dataset of known (i.e., experimentally determined) 3D structures
```{python}
# load packages

import numpy as np
import argparse
import math
import pandas as pd
import importlib
# import py file with common functions
import part1 as part1
part1.load()
```
Read pdb file as python dictionary 
calculate the distances between atoms
```{python}
importlib.reload(part1)
atoms = part1.read_pdb_file('./data/pdb/4gxy.pdb')
distances = part1.calculate_ca_distances(atoms)
```
calculate the observed frequencies into python dataframe
```{python}
observed_frequencies = part1.calculate_observed_frequencies(distances)
observed_frequencies
```

```{python}
importlib.reload(part1)
reference_frequencies = part1.calculate_frequencies(distances)
reference_frequencies
```
# Part 2
## plot the scoring profiles
```{python}
import part2 as part2
score_path = './data/scores/'
output_path = './data/figs/'
```
```{python}
part2.plot_distributions(score_path, output_path)
```
# Part 3
## Make predictions
```{python}
import part3 as part3
import part1 as part1
```
```{python}
pdb_file = './data/pdb/4gxy.pdb'
atoms = part1.read_pdb_file(pdb_file)
distances = part1.calculate_ca_distances(atoms)
gibbs_energy = part3.calculate_gibbs_free_energy(distances)
print("==="*10,"\nEstimated Gibbs Free Energy:", gibbs_energy,'\n',"==="*10)
```
