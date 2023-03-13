import pandas as pdimport numpy as npimport rdkitfrom rdkit import Chemfrom rdkit import DataStructs# Opening the reaction databasedata = pd.read_csv(r'/Users/maxbutterworth/Library/CloudStorage/OneDrive-TheUniversityofNottingham/Nottingham/Year 3/Synoptic Module/Research Project/Solvent-Selection/Test Database.csv')data = data.to_numpy() # Assigning Pandas database to numpy array# Defining category arraysreactants = data[0:, 0]solvents = data[0:, 1]smiles = data[0:, 6]ms = [] # Defining the molecular SMILES listfps = [] # Defining the fingerprints listsim_scores = [] # Defining a similarity scores list  # Defining a function to search the reactiondef reaction_search(input_smiles):    # Setting up the SMILES fingerprints in dataset    for i in range(0, len(smiles)):        #print(smiles[i])        ms.append(Chem.rdmolfiles.MolFromSmiles(smiles[i]))        fps.append(Chem.RDKFingerprint(ms[i]))    # Converting input smiles into a fingerprint    smi2mol = Chem.rdmolfiles.MolFromSmiles(input_smiles)    mol2fing = Chem.RDKFingerprint(smi2mol)        # Determining similarity scores    for j in range(0, len(fps)):        sim_scores.append(DataStructs.FingerprintSimilarity(mol2fing, fps[j]))        # Return the highest similarity score and the corresponding reaction SMILES    return np.max(sim_scores), smiles[np.argmax(np.max(sim_scores))]simscore, simreaction = reaction_search('COCCO')print(simscore)print(simreaction)def Solvent_Selection(reactant1, reactant2, product1, product2):        search_ms = [0,0,0,0]    search_fps = []        similarity_scores = []        # Setting up the fingerprints for input SMILES    search_ms[0] = Chem.rdmolfiles.MolFromSmiles(reactant1)    search_ms[1] = Chem.rdmolfiles.MolFromSmiles(reactant2)    search_ms[2] = Chem.rdmolfiles.MolFromSmiles(product1)    search_ms[3] = Chem.rdmolfiles.MolFromSmiles(product2)        for i in range(0,4):        search_fps.append(Chem.RDKFingerprint(search_ms[i]))            # Iterating over the dataset to compare similarity between dataset reactions and input SMILES    for j in range(0,len(reactants)):                sim_r1 = DataStructs.FingerprintSimilarity(search_fps[0], reactants1_fps[j])        sim_r2 = DataStructs.FingerprintSimilarity(search_fps[1], reactants2_fps[j])        sim_p1 = DataStructs.FingerprintSimilarity(search_fps[2], products1_fps[j])        sim_p2 = DataStructs.FingerprintSimilarity(search_fps[3], products2_fps[j])                similarity_scores[j] = sim_r1*sim_r2*sim_p1*sim_p2