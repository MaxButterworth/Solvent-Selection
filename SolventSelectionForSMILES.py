import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from ord_schema import message_helpers
from ord_schema.proto import dataset_pb2

# Opening the reaction database
data = message_helpers.load_message('ord_dataset-00005539a1e04c809a9a78647bea649c.pb.gz', dataset_pb2.Dataset)
df = message_helpers.messages_to_dataframe(data.reactions, drop_constant_columns=True) # Assigning data to a Pandas dataframe, df

# Defining an array for the reaction identifiers in ORD dataset
identify = df['identifiers[1].value']

# Defining reactant arrays for dataset. All elements are type str (see splicing a Pandas Data frame)
reactants1 = df['inputs["amine"].components[0].identifiers[0].value'] # Amine reactants
reactants2 = df['inputs["aryl halide"].components[0].identifiers[0].value'] # Aryl halide

# Defining product arrays for dataset. All elements are type str
products1 = df['outcomes[0].products[0].identifiers[0].value'] # Products

# Defining reagent arrays for dataset. All elements are type str
solvents = df['inputs["Solvent"].components[0].identifiers[0].value']
catalysts1 = df['inputs["metal and ligand"].components[0].identifiers[0].value'] # Assume this is a catalyst
catalysts2 = df['inputs["metal and ligand"].components[1].identifiers[0].value']
base = df['inputs["Base"].components[0].identifiers[0].value']
yields = df['outcomes[0].products[0].measurements[0].percentage.value']
temps = df['conditions.temperature.setpoint.value'] # Assume in degrees Celcius

# Defining empty arrays to append the fingerprints for the ORD dataset
reactants1_fps = []
reactants2_fps = []
products1_fps = []

# Defining empty arrays for sustainability metrics
sus_temp = [] # List to append sustainability ratings for reaction temperature
yieldband = [] # List to append yield sustainability bands

# Converting ORD dataset SMILES to fingerprints
for x in range(0, len(reactants1)):
    reactants1_ms = Chem.rdmolfiles.MolFromSmiles(reactants1[x])
    reactants1_fps.append( Chem.RDKFingerprint(reactants1_ms))
    
    reactants2_ms = Chem.rdmolfiles.MolFromSmiles(reactants2[x])
    reactants2_fps.append(Chem.RDKFingerprint(reactants2_ms))
    
    products1_ms = Chem.rdmolfiles.MolFromSmiles(products1[x])
    products1_fps.append(Chem.RDKFingerprint(products1_ms))
    
# Assessing sustainability based on temperature (using paper DOI: 10.1039/C5GC00340G)
for y in range(0, len(temps)):
    temp_value = temps[y]
    
    if np.isnan(temp_value)==True:
        sus_temp.append('No Data') # A value of zero means temperature data is not available
    
    elif 0 <= temp_value <= 70:
        sus_temp.append('Recommended') # A value of one is given a green label, i.e. most sustainable
    
    elif -20 <= temp_value < 0 or 70 < temp_value <= 140:
        sus_temp.append('Problematic') # A value of two is given a yellow label
    
    else:
        sus_temp.append('Hazardous') # A value of three is given a red label, i.e. least sustainable
        
for z in range(0, len(yields)):
    yield_value = yields[z]
    
    if yield_value > 0.89:
        yieldband.append('Green') # Given a green label, i.e. most sustainable
    elif 0.70 <= yield_value <= 0.89:
        yieldband.append('Amber') # Given a amber label
    elif yield_value < 0.70:
        yieldband.append('Red') # Given a red label, i.e. least sustainable
    else:
        yieldband.append('No Data') # Yield data is not available
  
# Defining a function to search the reaction
def Solvent_Selection(reactant1, reactant2, product1):
    
    # Setting up arrays
    search_ms = [0, 0, 0] # Mol from SMILES array
    search_fps = [] # Array for fingerprints to compare to database reaction fingerprints
    
    similarity_scores = [] # Array to append similarity scores
    k = 0 # Set up counter for top three similarity scores sorting loop
    top_5_sim = [] # Array for top three indices
    top_5_indices = [] # Array for indices of top three similarities
    
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
    
    # Finding the three highest similarity scores
    while k <=4:
        top_5_sim.append(np.max(similarity_scores))
        top_5_indices.append(np.argmax(similarity_scores))
        similarity_scores.remove(np.max(np.max(similarity_scores)))
        k += 1
        
    # Return the highest similarity score, the corresponding reaction SMILES and the correspoding solvent
    return top_5_sim, top_5_indices

def Predict_Reaction_Solvent(reactant1, reactant2, product1):
    sim, indices = Solvent_Selection(reactant1, reactant2, product1)
    
    x_labels = ['Reaction Similarity Score', 'Suggested Solvent', 'Temperature Sustainability Rating', 'Yield Sustainability Rating']
    y_labels = ['Reaction 1', 'Reaction 2', 'Reaction 3', 'Reaction 4', 'Reaction 5']
    
    results = [[round(sim[0], 5), solvents[indices[0]], sus_temp[indices[0]], yieldband[indices[0]]],
               [round(sim[1], 5), solvents[indices[1]], sus_temp[indices[1]], yieldband[indices[1]]],
               [round(sim[2], 5), solvents[indices[2]], sus_temp[indices[2]], yieldband[indices[2]]],
               [round(sim[3], 5), solvents[indices[3]], sus_temp[indices[3]], yieldband[indices[3]]],
               [round(sim[4], 5), solvents[indices[4]], sus_temp[indices[4]], yieldband[indices[4]]]]
    
    output = pd.DataFrame(results, y_labels, x_labels)
    
    output.style
    
    #print(output)
    
    return results


#Algorithm for a single reaction

# Defining a function to search the reaction
def Solvent_Selection_Single(reactant1, product1):
    
    # Setting up arrays
    search_ms = [0, 0] # Mol from SMILES array
    search_fps = [] # Array for fingerprints to compare to database reaction fingerprints
    
    similarity_scores = [] # Array to append similarity scores
    k = 0 # Set up counter for top three similarity scores sorting loop
    top_5_sim = [] # Array for top three indices
    top_5_indices = [] # Array for indices of top three similarities
    
    # Setting up the fingerprints for input SMILES
    search_ms[0] = Chem.rdmolfiles.MolFromSmiles(reactant1)
    search_ms[1] = Chem.rdmolfiles.MolFromSmiles(product1)
    
    for i in range(0,2):
        search_fps.append(Chem.RDKFingerprint(search_ms[i]))
        
    # Iterating over the dataset to compare similarity between dataset reactions and input SMILES
    for j in range(0,len(reactants1)):
        
        sim_r1 = DataStructs.FingerprintSimilarity(search_fps[0], reactants1_fps[j])
        sim_p1 = DataStructs.FingerprintSimilarity(search_fps[1], products1_fps[j])
        
        similarity_scores.append((sim_r1 + sim_p1)/2) # Linear average to determine overall similarity score
    
    # Finding the three highest similarity scores
    while k <=4:
        top_5_sim.append(np.max(similarity_scores))
        top_5_indices.append(np.argmax(similarity_scores))
        similarity_scores.remove(np.max(np.max(similarity_scores)))
        k += 1
        
    # Return the highest similarity score, the corresponding reaction SMILES and the correspoding solvent
    return top_5_sim, top_5_indices

def Predict_Reaction_Solvent_Single(reactant1, product1):
    sim, indices = Solvent_Selection_Single(reactant1, product1)
    
    x_labels = ['Reaction Similarity Score', 'Suggested Solvent', 'Temperature Sustainability Rating', 'Yield Sustainability Rating']
    y_labels = ['Reaction 1', 'Reaction 2', 'Reaction 3', 'Reaction 4', 'Reaction 5']
    
    results = [[round(sim[0], 5), solvents[indices[0]], sus_temp[indices[0]], yieldband[indices[0]]],
               [round(sim[1], 5), solvents[indices[1]], sus_temp[indices[1]], yieldband[indices[1]]],
               [round(sim[2], 5), solvents[indices[2]], sus_temp[indices[2]], yieldband[indices[2]]],
               [round(sim[3], 5), solvents[indices[3]], sus_temp[indices[3]], yieldband[indices[3]]],
               [round(sim[4], 5), solvents[indices[4]], sus_temp[indices[4]], yieldband[indices[4]]]]
    
    output = pd.DataFrame(results, y_labels, x_labels)
    
    output.style
    
    print(output)
    
    return results



# Allowing for user input of reaction SMILES and running the algorithm
#reactant1 = input('Input the first reactant in SMILES format: ')
#reactant2 = input('Input the second reactant in SMILES format: ')
#product1 = input('Input the product in SMILES format: ')

#Predict_Reaction_Solvent(reactant1, reactant2, product1)
#Predict_Reaction_Solvent_Single(reactant1, product1)

# Test input
#sim, indices = Solvent_Selection("CC(C)(C)C1=CC=C(C=C1)B(O)O", "BrC1=CC=CC=C1C=O", "CC(C)(C)C1=CC=C(C=C1)C1=C(C=O)C=CC=C1") # From this paper: https://doi.org/10.1021/acs.joc.1c00871





