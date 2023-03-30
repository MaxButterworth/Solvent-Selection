# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 09:53:04 2023

@author: pmyjb31
"""

yields = df['outcomes[0].products[0].measurements[0].percentage.value']

def Solvent_Selection(reactant1, reactant2, product1):
    
    # Setting up arrays
    search_ms = [0,0,0] # Mol from SMILES array
    search_fps = [] # Array for fingerprints to compare to database reaction fingerprints
    
    similarity_scores = [] # Array to append similarity scores
    k = 0 # Set up counter for top three similarity scores sorting loop
    top_3_sim = [] # Array for top three indices
    top_3_indices = [] # Array for indices of top three similarities
    yieldband = [] # Array for yield bands of all reactions
    top_3_yieldband = [] # Array for yield bands of top three similarities
    
    # Setting up the fingerprints for input SMILES
    search_ms[0] = Chem.rdmolfiles.MolFromSmiles(reactant1)
    search_ms[1] = Chem.rdmolfiles.MolFromSmiles(reactant2)
    search_ms[2] = Chem.rdmolfiles.MolFromSmiles(product1)
    #search_ms[3] = Chem.rdmolfiles.MolFromSmiles(product2)
    
    for i in range(0,3):
        search_fps.append(Chem.RDKFingerprint(search_ms[i]))
        
    # Iterating over the dataset to compare similarity between dataset reactions and input SMILES
    for j in range(0,len(reactants1)):
        
        sim_r1 = DataStructs.FingerprintSimilarity(search_fps[0], reactants1_fps[j])
        sim_r2 = DataStructs.FingerprintSimilarity(search_fps[1], reactants2_fps[j])
        sim_p1 = DataStructs.FingerprintSimilarity(search_fps[2], products1_fps[j])
        # sim_p2 = DataStructs.FingerprintSimilarity(search_fps[3], products2_fps[j])
        
        similarity_scores.append((sim_r1 + sim_r2 + sim_p1)/3) # Linear average to determine overall similarity score
    
        if yields[j] > 0.89:
            yieldband.append('green')
        elif yields[j] =< 0.89 and yields[j] >= 0.70:
            yieldband.append('amber')
        elif yields[j] < 0.70:
            yieldband.append('red')    
        else:
            yieldband.append('unknown')
        
        
    # Finding the three highest similarity scores
    while k <=2:
        top_3_sim.append(np.max(similarity_scores))
        top_3_indices.append(np.argmax(similarity_scores))
        similarity_scores.remove(np.max(np.max(similarity_scores)))
        top_3_yieldband.append(yieldband[top_3_indices[k]])
        k += 1
        
    # Return the highest similarity score, the corresponding reaction SMILES and the correspoding solvent
    return top_3_sim, top_3_indices,


#sim, indices = Solvent_Selection("CC(C)(C)C1=CC=C(C=C1)B(O)O", "BrC1=CC=CC=C1C=O", "CC(C)(C)C1=CC=C(C=C1)C1=C(C=O)C=CC=C1") # From this paper: https://doi.org/10.1021/acs.joc.1c00871

def Predict_Reaction_Solvent(reactant1, reactant2, product1):
    sim, indices, yieldband = Solvent_Selection(reactant1, reactant2, product1)
    print(f'The most similar reaction in the database has a similarity score of {round(sim[0], 5)}. The recommended solvent is {solvents[indices[0]]}. It has a yield banding of {yieldband[0]}.')
    
    print('')
    
    print(f'The second most similar reaction in the database has a similarity score of {round(sim[1], 5)}. The second recommended solvent is {solvents[indices[1]]}. It has a yield banding of {yieldband[1]}.')
    
    print('')
    
    print(f'The third most similar reaction in the database has a similarity score of {round(sim[2], 5)}. The third recommended solvent is {solvents[indices[2]]}. It has a yield banding of {yieldband[2]}.')   