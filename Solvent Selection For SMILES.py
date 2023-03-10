import pandas as pdimport numpy as npfrom rdkit import Chemfrom rdkit import DataStructsfrom ord_schema import message_helpersfrom ord_schema.proto import dataset_pb2# Opening the reaction databasedata = message_helpers.load_message('ord_dataset-00005539a1e04c809a9a78647bea649c.pb.gz', dataset_pb2.Dataset)df = message_helpers.messages_to_dataframe(data.reactions, drop_constant_columns=True) # Assigning data to a Pandas dataframe, df"""# Defining category arrays (see splicing a Pandas Data frame)reactants = df.iloc[]solvents = df.iloc[]smiles = df.iloc[]ms = [] # Defining the molecular SMILES listfps = [] # Defining the fingerprints listsim_scores = [] # Defining a similarity scores list  # Defining a function to search the reactiondef reaction_search(input_smiles):    # Setting up the SMILES fingerprints in dataset    for i in range(0, len(smiles)):        #print(smiles[i])        ms.append(Chem.rdmolfiles.MolFromSmiles(smiles[i]))        fps.append(Chem.RDKFingerprint(ms[i]))    # Converting input smiles into a fingerprint    smi2mol = Chem.rdmolfiles.MolFromSmiles(input_smiles)    mol2fing = Chem.RDKFingerprint(smi2mol)        # Determining similarity scores    for j in range(0, len(fps)):        sim_scores.append(DataStructs.FingerprintSimilarity(mol2fing, fps[j]))        # Return the highest similarity score and the corresponding reaction SMILES    return np.max(sim_scores), smiles[np.argmax(np.max(sim_scores))]simscore, simreaction = reaction_search('COCCO')print(simscore)print(simreaction)"""