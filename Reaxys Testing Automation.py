import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from ord_schema import message_helpers
from ord_schema.proto import dataset_pb2

import SolventSelectionForSMILES as SolSel # Import our solvent selection algorithm

REAXYS_reaction_SMILES = [] # Import reaction SMILES from REAXYS
test_results_sim = [] # Array for similarity of first reaction suggestions
test_results_solvent = [] # Array for solvent of first reaction suggestions
test_results_temp = [] # Array for temperature of first reaction suggestions
test_results_yield = [] # Array for yield of first reaction suggestions

for x in range(0, len(REAXYS_reaction_SMILES)):
    
    reaction_SMILES = REAXYS_reaction_SMILES[x]
    
    if "." not in reaction_SMILES: #Single reactant reaction
        for y in reaction_SMILES: # Looping through each character
            
            if y == ">": # y and y+1 are the reaction arrow
                reaction_arrow = y 
                break
        
        reactant1 = reaction_SMILES[0:reaction_arrow-1] # Reactant is from the first character up to the reaction arrow
        product1 = reaction_SMILES[reaction_arrow+2:len(reaction_SMILES)-1] # Product is from the reaction arrow up to the last character
        
        reaction_results = SolSel.Predict_Reaction_Solvent_Single(reactant1, product1)
        test_results_sim.append(reaction_results[0,0])
        test_results_solvent.append(reaction_results[0,1])
        test_results_temp.append(reaction_results[0,2])
        test_results_yield.append(reaction_results[0,3])
    
    else:
        for y in reaction_SMILES: # Looping through each character
            
            if y == ".": # y is the separator between reactant 1 and 2
                reactant_separator = y
            
            elif y == ">": # y and y+1 are the reaction arrow
                reaction_arrow = y
                break
        
        reactant1 = reaction_SMILES[0:reactant_separator-1] # reactant1 is from the first character up to the reactant_separator
        reactant2 = reaction_SMILES[reactant_separator+1:reaction_arrow-1]# reactant2 is from the reactant_separator up to the reaction_arrow
        product1 = reaction_SMILES[reaction_arrow+2:len(reaction_SMILES)-1] # product1 is from the reaction arrow up to the last character
        
        reaction_results = SolSel.Predict_Reaction_Solvent(reactant1, reactant2, product1)
        test_results_sim.append(reaction_results[0,0])
        test_results_solvent.append(reaction_results[0,1])
        test_results_temp.append(reaction_results[0,2])
        test_results_yield.append(reaction_results[0,3])
        
        
        
        

