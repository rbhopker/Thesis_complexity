#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 15:36:13 2022

@author: ricardobortothopker
"""

from pysmiles import read_smiles
import networkx as nx
import matplotlib.pyplot as plt
from complexFuncs import complexity
import numpy as np
    
# smiles = 'C12=C3C4=C5C6=C1C7=C8C9=C1C%10=C%11C(=C29)C3=C2C3=C4C4=C5C5=C9C6=C7C6=C7C8=C1C1=C8C%10=C%10C%11=C2C2=C3C3=C4C4=C5C5=C%11C%12=C(C6=C95)C7=C1C1=C%12C5=C%11C4=C3C3=C5C(=C81)C%10=C23'
Strychnine = "C1CN2CC3=CCOC4CC(=O)N5C6C4C3CC2C61C7=CC=CC=C75"
PubChemCS = "CC1=C(C=C(C=C1)Cl)N2CCN(CC2)C(=O)C3CC4=C(C=CC(=C4)OC)OC3"
Glucose = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"
Anandamide = "CCCCC/C=C\C/C=C\C/C=C\C/C=C\CCCC(=O)NCCO"
Adamantane = "C1C2CC3CC1CC(C2)C3"

Strychnine_mol,PubChemCS_mol,Glucose_mol = read_smiles(Strychnine),read_smiles(PubChemCS),read_smiles(Glucose)
Anandamide_mol,Adamantane_mol = read_smiles(Anandamide),read_smiles(Adamantane)

mols = [Strychnine_mol,PubChemCS_mol,Glucose_mol,Anandamide_mol,Adamantane_mol]
A = []
elements =[]
C = []
for mol in mols:
    fig,ax = plt.subplots()
    A.append(nx.to_numpy_matrix(mol))
    np.fill_diagonal(A[-1],1)
    elements.append(nx.get_node_attributes(mol, name = "element"))
    # nx.draw(mol, with_labels=True, labels = elements[-1], pos=nx.spring_layout(mol))
    # plt.gca().set_aspect('equal')
    C.append(complexity(A[-1]))