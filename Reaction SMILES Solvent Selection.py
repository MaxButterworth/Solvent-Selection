# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 01:45:19 2023

@author: joela
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit import DataStructs

def Solvent_Selection(reactant1, reactant2, product1, product2):
    
    search_ms = [0,0,0,0]
    search_fps = []
    
    similarity_scores = []
    
    search_ms[0] = Chem.rdmolfiles.MolFromSmiles(reactant1)
    search_ms[1] = Chem.rdmolfiles.MolFromSmiles(reactant2)
    search_ms[2] = Chem.rdmolfiles.MolFromSmiles(product1)
    search_ms[3] = Chem.rdmolfiles.MolFromSmiles(product2)
    
    for i in range(0,4):
        search_fps.append(Chem.RDKFingerprint(search_ms[i]))
    
    for j in range(0,len(reactants)):
        
        sim_r1 = DataStructs.FingerprintSimilarity(search_fps[0], reactants1_fps[j])
        sim_r2 = DataStructs.FingerprintSimilarity(search_fps[1], reactants2_fps[j])
        sim_p1 = DataStructs.FingerprintSimilarity(search_fps[2], products1_fps[j])
        sim_p2 = DataStructs.FingerprintSimilarity(search_fps[3], products2_fps[j])
        
        similarity_scores[j] = sim_r1*sim_r2*sim_p1*sim_p2