import pandas as pdimport numpy as npfrom rdkit import Chemfrom rdkit import DataStructs# Opening the reaction databasedata = pd.read_csv(r'/Users/maxbutterworth/Library/CloudStorage/OneDrive-TheUniversityofNottingham/Nottingham/Year 3/Synoptic Module/Research Project/Solvent-Selection/Database copy.csv')data = data.to_numpy() # Assigning Pandas database to numpy array# Defining category arraysreactants = data[0:, 0]solvents = data[0:, 1]smiles = data[0:, 6]ms = [] # Defining the molecular SMILES listfps = [] # Defining the fingerprints listsim_scores = [] # Defining a similarity scores list  # Defining a function to search the reactiondef reaction_search(input_smiles):    # Setting up the SMILES fingerprints in dataset    for i in range(0, len(smiles)):        #print(smiles[i])        ms.append(Chem.rdmolfiles.MolFromSmiles(smiles[i]))        fps.append(Chem.RDKFingerprint(ms[i]))    # Converting input smiles into a fingerprint    smi2mol = Chem.rdmolfiles.MolFromSmiles(input_smiles)    mol2fing = Chem.RDKFingerprint(smi2mol)        # Determining similarity scores    for j in range(0, len(fps)):        sim_scores.append(DataStructs.FingerprintSimilarity(mol2fing,fps[j]))        # Return the highest similarity score and the corresponding reaction SMILES    return np.max(sim_scores), smiles[np.argmax(np.max(sim_scores))]simscore, simreaction = reaction_search('COCCO')print(simscore)print(simreaction)