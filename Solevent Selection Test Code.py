import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit import DataStructs
import plotly as px
import ord_schema



# Opening the reaction database
database = pd.read_csv(r'/Users/maxbutterworth/Library/CloudStorage/OneDrive-TheUniversityofNottingham/Nottingham/Year 3/Synoptic Module/Research Project/Solvent-Selection/Test Database.csv')
database = database.to_numpy() # Assigning Pandas database to numpy array

# Defining category arrays
reagents = database[0:, 0]
solvents = database[0:, 1]
smiles = database[0:, 6]
ms = [] # Defining the molecular SMILES list
fps = [] # Defining the fingerprints list
#print(smiles)

# Setting up the SMILES fingerprints
for i in range(0, len(smiles)):
    mol_from_smiles = Chem.rdmolfiles.MolFromSmiles((f'{smiles[i]}'))
    #ms.append(f'{mol_from_smiles}')
    fingerprint = Chem.RDKFingerprint(f'{mol_from_smiles}')
    fps.append(f'{fingerprint}')
    #print(type(ms))

#print(fps)    
# Defining a function to search the reaction
#def reaction_search(input_smiles):
    
    

#print(database)