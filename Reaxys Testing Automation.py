import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from ord_schema import message_helpers
from ord_schema.proto import dataset_pb2

import SolventSelectionForSMILES as SolSel # Import our solvent selection algorithm

df = pd.read_excel('Reaxys Test Data - Copy.xlsx')

REAXYS_reaction_SMILES = df['Reaction'].values # Import reaction SMILES from REAXYS

test_results_sim = [] # Array for similarity of first reaction suggestions
test_results_solvent = [] # Array for solvent of first reaction suggestions
test_results_temp = [] # Array for temperature of first reaction suggestions
test_results_yield = [] # Array for yield of first reaction suggestions

for x in range(0, len(REAXYS_reaction_SMILES)):
    
    reaction_SMILES = REAXYS_reaction_SMILES[x]
    
    y_number = -1
    
    if "." not in reaction_SMILES: #Single reactant reaction
        for y in reaction_SMILES: # Looping through each character
            y_number += 1 
        
            if y == ">": # y and y+1 are the reaction arrow
                reaction_arrow = y_number 
                break
        
        reactant1 = reaction_SMILES[0:reaction_arrow] # Reactant is from the first character up to the reaction arrow
        product1 = reaction_SMILES[reaction_arrow+2:len(reaction_SMILES)] # Product is from the reaction arrow up to the last character
        
        reaction_results = SolSel.Predict_Reaction_Solvent_Single(reactant1, product1)
        
        test_results_sim.append(reaction_results[0][0])
        test_results_solvent.append(reaction_results[0][1])
        test_results_temp.append(reaction_results[0][2])
        test_results_yield.append(reaction_results[0][3])
        
    
    else:
        for y in reaction_SMILES: # Looping through each character
            
            y_number += 1
            
            if y == ".": # y is the separator between reactant 1 and 2
                reactant_separator = y_number
            
            elif y == ">": # y and y+1 are the reaction arrow
                reaction_arrow = y_number
                break
        
        reactant1 = reaction_SMILES[0:reactant_separator] # reactant1 is from the first character up to the reactant_separator
        reactant2 = reaction_SMILES[reactant_separator+1:reaction_arrow]# reactant2 is from the reactant_separator up to the reaction_arrow
        product1 = reaction_SMILES[reaction_arrow+2:len(reaction_SMILES)] # product1 is from the reaction arrow up to the last character
        
        reaction_results = SolSel.Predict_Reaction_Solvent(reactant1, reactant2, product1)
        test_results_sim.append(reaction_results[0][0])
        test_results_solvent.append(reaction_results[0][1])
        test_results_temp.append(reaction_results[0][2])
        test_results_yield.append(reaction_results[0][3])
        
    #print (test_results_sim[x])
    
from openpyxl import Workbook

def export_test_data(path):
    workbook = Workbook()
    sheet = workbook.active
    
    sheet["A1"] = "Reaction Similarity"
    sheet["B1"] = "Suggested Solvent SMILES"
    sheet["C1"] = "Reaction Temperature Sustainability"
    sheet["D1"] = "Reaction Yield Sustainability"
    
    for i in range(2, len(test_results_sim)+2):
        sheet["A"+str(i)] = str(test_results_sim[i-2])
        sheet["B"+str(i)] = str(test_results_solvent[i-2])
        sheet["C"+str(i)] = str(test_results_temp[i-2])
        sheet["D"+str(i)] = str(test_results_yield[i-2])
    
    workbook.save(path)
        
export_test_data("Test Results.xlsx")        
        
        

